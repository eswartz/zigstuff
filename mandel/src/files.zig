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

        try self.load_version_1_1(reader);
    }

    fn load_version_1_1(self: *Self, reader: File.Reader) !void {
        var cxStr = try self.read_string(reader);
        var cyStr = try self.read_string(reader);

        self.params.zoom = try std.fmt.parseInt(u16, try self.read_string(reader), 10);

        // self.storage.alloc.free(self.params.cx);
        // self.storage.alloc.free(self.params.cy);

        var bits = self.params.getDefaultIntSize() << 6;
        self.params.words = bits >> 6;
        self.params.cx = try self.parse_decimal(cxStr, bits);
        self.params.cy = try self.parse_decimal(cyStr, bits);

        var minIters = try std.fmt.parseInt(u32, try self.read_string(reader), 10);
        _ = minIters;
        self.params.iters = try std.fmt.parseInt(u32, try self.read_string(reader), 10);

        self.params.sx = try std.fmt.parseInt(u32, try self.read_string(reader), 10);
        self.params.sy = try std.fmt.parseInt(u32, try self.read_string(reader), 10);

        // const ps = self.params.blockSize;
        // const SEGS: u32 = @divExact(std.math.min(params.sx, params.sy), ps);
        // const floatStep = BigFloatType.fromFloat(params.span() / @intToFloat(bignum.NativeFloat, ps * SEGS));

        var layer = try self.storage.ensure(self.params.zoom, self.params.blockSize);

        var zstr = try std.compress.gzip.gzipStream(self.arena, reader);
        var exp = self.params.sx * self.params.sy;
        const buf = try zstr.reader().readAllAlloc(self.arena, exp*4);
        const iters = @ptrCast([*]align(1) i32, buf)[0..exp];

        const blockSize = self.params.blockSize;
        const floatStep = self.params.pixelSpan();
        var px : u32 = 0;
        var py : u32 = 0;
        var ix : i1024 = 0;
        var iy : i1024 = 0;
        var block : *mandel.BlockData = undefined;
        for (iters) |it| {
            // std.debug.print("{x}{s}", .{ @byteSwap(it), if (c%16 == 15) "\n" else " " }); c+= 1;
            if ((px & (blockSize-1)) == 0) {
                block = try layer.ensure(mandel.BigCoord.init(ix, iy));
            }
            block.setIter(px & (blockSize-1), py & (blockSize-1), @byteSwap(it));
            px += 1;
            ix += floatStep;
            if (px >= self.params.sx) {
                px = 0;
                py += 1;
                iy += floatStep;
            }
        }
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
        return l;
    }

    /// Read a normal long decimal and convert to our hex format (.XXXX.XXX.XXX)
    fn parse_decimal(self: Self, str: []const u8, bits: u16) ![]u8 {
        var mpf : gmp.mpf_t = undefined;
        gmp.mpf_init(&mpf);
        defer gmp.mpf_clear(&mpf);

        gmp.mpf_set_prec(&mpf, bits + 8);
        if (gmp.mpf_set_str(&mpf, str.ptr, 10) != 0)
            return error.InvalidFloat;

        // scale to fixed int size, where the top 8 bits are the integral part
        gmp.mpf_mul_2exp(&mpf, &mpf, bits - 8);
        {
            var dec : gmp.mp_exp_t = undefined;
            var strptr = gmp.mpf_get_str(null, &dec, 16, 256, &mpf);
            std.debug.print("xx> {s}.{s}\n", .{strptr[0..@intCast(usize, dec)], strptr[@intCast(usize, dec)+1..strlen(strptr)]});
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

        var bits = (self.params.getDefaultIntSize() << 6);
        const cxStr = try self.write_decimal(self.params.cx, bits);
        const cyStr = try self.write_decimal(self.params.cy, bits);
        try self.write_string(writer, cxStr);
        try self.write_string(writer, cyStr);

        try self.write_string(writer, try std.fmt.allocPrint(self.arena, "{}", .{ self.params.zoom }));

        try self.write_string(writer, try std.fmt.allocPrint(self.arena, "{}", .{ self.params.iters }));
        try self.write_string(writer, try std.fmt.allocPrint(self.arena, "{}", .{ self.params.iters }));

        try self.write_string(writer, try std.fmt.allocPrint(self.arena, "{}", .{ self.params.sx }));
        try self.write_string(writer, try std.fmt.allocPrint(self.arena, "{}", .{ self.params.sy }));

        // // expect compression to be smaller than this...
        // var zcontent = try reader.readAllAlloc(self.arena, self.params.sx * self.params.sy * 4);
        // std.debug.print("zcontent size={}\n", .{zcontent.len});
    // _iters = content.decompress(_psx * _psy * 4, File.COMPRESSION_GZIP)
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
