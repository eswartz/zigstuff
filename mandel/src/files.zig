const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;

const bignum = @import("bignum.zig");
const Mandel = bignum.Mandel;

const mandel = @import("mandel.zig");
const Params = mandel.Params;

const File = std.fs.File;

const cimports = @import("cimports.zig");
const gmp = cimports.gmp;
const zlib = cimports.zlib;

var arena_instance = std.heap.ArenaAllocator.init(std.heap.raw_c_allocator);

pub const RenderedFile = struct {
    const Self = @This();

    params: *Params,
    storage : *mandel.MandelStorage,
    arena: Allocator,

    pub fn init(params: *Params, storage: *mandel.MandelStorage) !Self {
        return Self{
            .params = params,
            .storage = storage,
            .arena = arena_instance.allocator()
        };
    }

    pub fn deinit(self: *Self) void {
        _ = self;
    }

    pub fn load(self: *Self, path: []const u8) !void {
        var file = try std.fs.cwd().openFile(path, .{});
        defer file.close();

        var reader = file.reader();
        var version = try self.read_string(reader);
        if (!std.mem.eql(u8, version, "jmandel_1.1")) {
            return error.UnsupportedVersion;
        }

        try self.load_version_1_1(reader, (try file.stat()).size);
    }

    inline fn decompress_blocks(self : *Self, comptime BigIntType : type, psx : u32, psy : u32, iters: []const i32) !void {
        _ = psy;
        var workLoad = try Params.BlockWorkMaker(BigIntType).init(self.storage.alloc, self.params.*, self.storage);
        defer workLoad.deinit(self.storage.alloc);

        std.debug.print("iters={}, psx={}, psy={}\n", .{ iters.len, self.params.sx, self.params.sy });

        const ps = self.params.blockSize;
        var c : u32 = 0;
        for (workLoad.blocks) |block| {
            const px = block.px;
            const py = block.py;
            // std.debug.print("block={*}\n", .{ block.data });
            var oy : u32 = 0;
            while (oy < ps) : (oy += 1) {
                var ox : u32 = 0;
                while (ox < ps) : (ox += 1) {
                    const iterBE : i32 = iters[(py + oy) * psx + (px + ox)];
                    const iter = @byteSwap(iterBE);
                    block.data.setIter(ox, oy, iter);
                    c += 1;
                }
            }
        }
    }

    fn load_version_1_1(self: *Self, reader: File.Reader, fileSize: u64) !void {
        var cxStr = try self.read_string(reader);
        var cyStr = try self.read_string(reader);

        self.params.zoom = try std.fmt.parseInt(u16, try self.read_string(reader), 10);

        var bits = self.params.getDefaultIntSize() << 6;
        self.params.words = bits >> 6;
        self.params.cx = try self.parse_decimal(cxStr, bits);
        self.params.cy = try self.parse_decimal(cyStr, bits);

        var minIters = try std.fmt.parseInt(u32, try self.read_string(reader), 10);
        _ = minIters;
        self.params.iters = try std.fmt.parseInt(u32, try self.read_string(reader), 10);

        const psx = try std.fmt.parseInt(u32, try self.read_string(reader), 10);
        const psy = try std.fmt.parseInt(u32, try self.read_string(reader), 10);
        self.params.sx = psx;
        self.params.sy = psy;

        const exp : u32 = psx * psy;
        var iters = try self.arena.alloc(i32, exp);
        std.mem.set(i32, iters, 0);
        try self.decompress_iters(iters, reader, fileSize);

        var i : i32 = 0; var l : i32 = 0;
        for (iters) |iter| {
            if (iter != l) {
                std.debug.print("{} ", .{iter});
                l = iter;
                i += 1;
                if (@rem(i, 16) == 0) std.debug.print("\n", .{});
            }
        }

        try switch (self.params.words) {
            inline 1...bignum.MAXWORDS => |w| self.decompress_blocks(bignum.BigInt(w << 6), psx, psy, iters),
            else => unreachable,
        };
    }

    fn decompress_iters(self: Self, iters : []i32, reader : File.Reader, fileSize: u64) !void {
        const zmem = try self.arena.alloc(u8, fileSize);
        const zlen = try reader.readAll(zmem);

        var zstr : zlib.z_stream = undefined;
        std.mem.set(u8, @ptrCast([*]u8, &zstr)[0..@sizeOf(zlib.z_stream)], 0);
        // window size (15, default) gzip header (+32)
        if (zlib.inflateInit2(&zstr, 15 + 32) != zlib.Z_OK) return error.ZLibError;
        zstr.next_in = zmem.ptr;
        zstr.avail_in = @intCast(c_uint, zlen);

        var i8ptr : usize = 0;
        var iters8 = @ptrCast([*]u8, iters);

        while (i8ptr < zlen) {
            zstr.avail_out = @intCast(c_uint, fileSize);
            zstr.next_out = iters8 + i8ptr;

            const res = zlib.inflate(&zstr, if (zstr.avail_in == 0) zlib.Z_FINISH else zlib.Z_NO_FLUSH);
            if (res == zlib.Z_STREAM_END) break;
            if (res == zlib.Z_DATA_ERROR) {
                std.debug.print("zlib error: {s}\n", .{zstr.msg});
                return error.FileFormat;
            }

            i8ptr += zstr.total_out;
            zstr.total_out = 0;
        }

        _ = zlib.inflateEnd(&zstr);
    }

    fn read_string(self: Self, reader : File.Reader) ![]u8 {
        var length = try reader.readIntBig(u16);
        var str = try self.arena.alloc(u8, length);
        var len = try reader.read(str);
        if (len != length) return error.BadFileFormat;
        std.debug.print("--> {s}\n", .{str});
        return str;
    }

    fn strlen(s:[*c]u8) usize {
        var l : usize = 0;
        while (s[l] != 0) l += 1;
        std.debug.print("len: {}\n", .{l});
        return l;
    }

    /// Read a normal long decimal and convert to our hex format (.XXXX.XXX.XXX)
    fn parse_decimal(self: Self, str: []const u8, bits: u16) ![]u8 {
        var mpf : gmp.mpf_t = undefined;
        gmp.mpf_init(&mpf);
        defer gmp.mpf_clear(&mpf);

        gmp.mpf_set_prec(&mpf, bits + 8);
        _ = gmp.mpf_set_str(&mpf, str.ptr, 10);

        // scale to fixed int size, where the top 8 bits are the integral part
        gmp.mpf_mul_2exp(&mpf, &mpf, bits - 8);
        {
            var exp : gmp.mp_exp_t = undefined;
            const strptr = gmp.mpf_get_str(null, &exp, 16, 0, &mpf);
            const slen = if (exp != 0 and strptr != null) strlen(strptr) else 0;
            if (exp == 0)
                std.debug.print("xx> 0\n", .{})
            else if (exp < slen)
                std.debug.print("xx> {s}.{s}\n", .{strptr[0..@intCast(usize, exp)], strptr[@intCast(usize, exp)+1..slen]})
            else
                std.debug.print("xx> {s}.0+\n", .{strptr[0..slen]});
        }

        // make into int
        var mpz : gmp.mpz_t = undefined;
        gmp.mpz_init(&mpz);
        defer gmp.mpz_clear(&mpz);

        gmp.mpz_set_f(&mpz, &mpf);

        // now get data
        var strptr = gmp.mpz_get_str(null, 16, &mpz);
        const istr = strptr[0..strlen(strptr)];
        std.debug.print("##> {s}\n", .{istr});

        var list = std.ArrayList(u8).init(self.storage.alloc);
        var p : usize = istr.len;
        if (p % 16 != 0) {
            // partial first word, pad with zeroes
            var q : isize = @mod(-@intCast(isize, p), 16);
            try list.append('.');
            while (q != 0) : (q -= 1) {
                try list.append('0');
            }
            p %= 16;
            try std.fmt.format(list.writer(), "{s}", .{istr[0..p]});
        }
        // full words
        while (p + 16 < istr.len) {
            try std.fmt.format(list.writer(), ".{s}", .{istr[p..p+16]});
            p += 16;
        }
        if (p < istr.len) {
            // partial last word
            var q = istr.len - p;
            try std.fmt.format(list.writer(), ".{s}", .{istr[p..p + q]});
            while (q < 16) : (q -= 1) {
                try list.append('0');
            }
        }

        const bistr = list.toOwnedSlice();
        std.debug.print("ii> {s}\n", .{bistr});

        return bistr;
    }

    pub fn save(self: *Self, path: []const u8) !void {
        var file = try std.fs.cwd().createFile(path, .{.truncate = true});
        defer file.close();

        var writer = file.writer();
        try self.write_string(writer, "jmandel_1.1");

        // var bits = self.params.words << 6;
        var bits = self.params.getDefaultIntSize() << 6;

        const cxStr = try self.write_decimal(self.params.cx, bits);
        const cyStr = try self.write_decimal(self.params.cy, bits);
        try self.write_string(writer, cxStr);
        try self.write_string(writer, cyStr);

        try self.write_string(writer, try std.fmt.allocPrint(self.arena, "{}", .{ self.params.zoom }));

        try self.write_string(writer, try std.fmt.allocPrint(self.arena, "{}", .{ self.params.iters }));
        try self.write_string(writer, try std.fmt.allocPrint(self.arena, "{}", .{ self.params.iters }));

        try self.write_string(writer, try std.fmt.allocPrint(self.arena, "{}", .{ self.params.sx }));
        try self.write_string(writer, try std.fmt.allocPrint(self.arena, "{}", .{ self.params.sy }));

        const iters = try self.arena.alloc(i32, self.params.sx * self.params.sy);
        std.mem.set(i32, iters, 0);

        try switch (self.params.words) {
            inline 1...bignum.MAXWORDS => |w| self.compress_blocks(bignum.BigInt(w << 6), self.params.sx, self.params.sy, iters),
            else => unreachable,
        };

        try self.compress_iters(iters, writer);

        // var zstr = try std.compress.gzip.gzipStream(self.arena, reader);
        // var exp = psx * psy;

        // var iters : []i32 = try self.arena.alloc(i32, exp);
        // const end = try zstr.reader().readAll(@ptrCast([*]u8, iters)[0..exp*4]);
        // if (end < exp * 4) return error.FileFormatError;

        // try switch (self.params.words) {
        //     inline 1...bignum.MAXWORDS => |w| self.decompress_blocks(bignum.BigInt(w << 6), psx, psy, iters),
        //     else => unreachable,
        // };

        // // expect compression to be smaller than this...
        // var zcontent = try reader.readAllAlloc(self.arena, self.params.sx * self.params.sy * 4);
        // std.debug.print("zcontent size={}\n", .{zcontent.len});
    // _iters = content.decompress(_psx * _psy * 4, File.COMPRESSION_GZIP)
    }

    fn compress_blocks(self : *Self, comptime BigIntType : type, psx : u32, psy : u32, iters: []i32) !void {
        _ = psy;
        var workLoad = try Params.BlockWorkMaker(BigIntType).init(self.storage.alloc, self.params.*, self.storage);
        defer workLoad.deinit(self.storage.alloc);

        std.debug.print("iters={}, psx={}, psy={}\n", .{ iters.len, self.params.sx, self.params.sy });

        const ps = self.params.blockSize;
        var c : u32 = 0;
        for (workLoad.blocks) |block| {
            const px = block.px;
            const py = block.py;
            // std.debug.print("block={*}\n", .{ block.data });
            var oy : u32 = 0;
            while (oy < ps) : (oy += 1) {
                var ox : u32 = 0;
                while (ox < ps) : (ox += 1) {
                    const iter : i32 = block.data.iter(ox, oy);
                    const iterBE = @byteSwap(iter);
                    iters[(py + oy) * psx + (px + ox)] = iterBE;
                    c += 1;
                }
            }
        }
    }

    fn compress_iters(self: Self, iters: []i32, writer: File.Writer) !void {
        var obuf = try self.arena.alloc(u8, 4096);

        var zstr : zlib.z_stream = undefined;
        std.mem.set(u8, @ptrCast([*]u8, &zstr)[0..@sizeOf(zlib.z_stream)], 0);

        // var i8ptr : usize = 0;
        var i8len : usize = iters.len * 4;

        // zstr.next_in = @ptrCast([*c] u8, @alignCast(1, iters.ptr));
        zstr.next_in = @ptrCast([*] u8, iters);
        zstr.avail_in = @intCast(c_uint, i8len);

        // window size (15, default) gzip header (+16)
        const ires = zlib.deflateInit2(&zstr, zlib.Z_DEFAULT_COMPRESSION, zlib.Z_DEFLATED, 15 + 16, 8, zlib.Z_DEFAULT_STRATEGY);
        if (ires != zlib.Z_OK) {
            std.debug.print("init error: {}\n", .{ires});
            return error.ZLibError;
        }

        while (zstr.avail_in != 0) {
            zstr.avail_out = @intCast(c_uint, obuf.len);
            zstr.next_out = obuf.ptr; //@ptrCast([*c] u8, &obuf);

            const res = zlib.deflate(&zstr, if (zstr.avail_in == 0) zlib.Z_FINISH else zlib.Z_NO_FLUSH);
            if (res == zlib.Z_STREAM_END) break;
            if (res == zlib.Z_DATA_ERROR) {
                std.debug.print("zlib error: {s}\n", .{zstr.msg});
                return error.FileFormat;
            }

            const olen = obuf.len - zstr.avail_out;
            std.debug.print("olen={}\n", .{olen});
            _ = try writer.write(obuf[0..olen]);
            // zstr.next_in += zstr.total_in;
            // zstr.total_out = 0;
        }

        _ = zlib.deflateEnd(&zstr);
    }

    fn write_string(self: Self, writer : File.Writer, str: []const u8) !void {
        _=self;
        try writer.writeIntBig(u16, @intCast(u16, str.len));
        try writer.writeAll(str);
    }

    /// Write our hex format (.XXXX.XXX.XXX) to "normal" decimal
    fn write_decimal(self: Self, istr: []const u8, bits: u16) ![]u8 {
        var mpz : gmp.mpz_t = undefined;
        gmp.mpz_init(&mpz);
        defer gmp.mpz_clear(&mpz);

        var p : usize = 0;
        while (p < istr.len) : (p += 17) {
            var w : u64 = try std.fmt.parseInt(u64, istr[p+1..p+17], 16);
            gmp.mpz_mul_2exp(&mpz, &mpz, 64);
            gmp.mpz_add_ui(&mpz, &mpz, w);
        }
        std.debug.print("{} - {}\n", .{p, istr.len});
        if (p != istr.len) unreachable;

        var mpf : gmp.mpf_t = undefined;
        gmp.mpf_init(&mpf);
        defer gmp.mpf_clear(&mpf);
        gmp.mpf_set_prec(&mpf, bits + 8);

        // get the int into the float
        gmp.mpf_set_z(&mpf, &mpz);

        // divide down
        gmp.mpf_div_2exp(&mpf, &mpf, bits - 8);

        var dec : gmp.mp_exp_t = undefined;
        var strptr = gmp.mpf_get_str(null, &dec, 10, 256, &mpf);
        std.debug.print("xx> ({}){s}\n", .{dec, strptr[0..strlen(strptr)]});

        var list = std.ArrayList(u8).init(self.arena);
        var i : usize = 0;
        if (dec < 0) {
            try list.append('.');
            while (dec < 0) : (dec += 1) {
                try list.append('0');
            }
        } else {
            try list.appendSlice(strptr[0..@intCast(usize, i)]);
            i = 0;
            try list.append('.');
            if (dec == 0) {
                try list.append('0');
            }
        }
        try list.appendSlice(strptr[i..strlen(strptr)]);

        const bistr = list.toOwnedSlice();
        std.debug.print("dd> {s}\n", .{bistr});
        return bistr;
    }
};

/////////////////

fn run_test(comptime callable: anytype) !void {
    const alloc = std.testing.allocator;
    var storage = try mandel.MandelStorage.init(alloc);
    defer storage.deinit();
    var params = try mandel.Params.init(alloc);
    defer params.deinit(alloc);

    try callable.doit(&storage, &params);

}

test "load1" {
    try run_test(struct{
        pub fn doit(storage: *mandel.MandelStorage, params: *mandel.Params) !void {
            var file = try RenderedFile.init(params, storage);
            defer file.deinit();
            try file.load("data/zoom10_000.dat");
        }
    });
}
test "save1" {
    try run_test(struct{
        pub fn doit(storage: *mandel.MandelStorage, params: *mandel.Params) !void {
            var file = try RenderedFile.init(params, storage);
            defer file.deinit();
            try file.load("data/zoom10_400.dat");
            try file.save("data/zoom10_400.dat.bak");
        }
    });
}
