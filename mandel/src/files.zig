const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;

const cimports = @import("cimports.zig");
const gmp = cimports.gmp;
const zlib = cimports.zlib;

const config = @import("config.zig");

const bignum = @import("bignum.zig");
const Mandel = bignum.Mandel;

const algo = @import("algo.zig");
const MandelParams = algo.Params;
const cache = @import("cache.zig");
const MandelStorage = cache.MandelStorage;

const File = std.fs.File;

pub const DecimalMath = struct {
    const Self = @This();

    fn strlen(s:[*c]u8) usize {
        var l : usize = 0;
        while (s[l] != 0) l += 1;
        std.debug.print("len: {}\n", .{l});
        return l;
    }

    /// Read a normal long decimal and convert to our hex format (.XXXX.XXX.XXX)
    /// caller owns returned array
    pub fn parse_decimal(alloc: Allocator, str: [:0]const u8, bits: u16) ![]u8 {
        var mpf : gmp.mpf_t = undefined;
        gmp.mpf_init(&mpf);
        defer gmp.mpf_clear(&mpf);

        // gmp.mpf_set_prec(&mpf, bits + 8);
        gmp.mpf_set_prec(&mpf, str.len * 4 + 8);
        _ = gmp.mpf_set_str(&mpf, str.ptr, 10);

        // scale to fixed int size, where the top 8 bits are the integral part
        gmp.mpf_mul_2exp(&mpf, &mpf, bits - 8);

        // make into int
        var mpz : gmp.mpz_t = undefined;
        gmp.mpz_init(&mpz);
        defer gmp.mpz_clear(&mpz);

        gmp.mpz_set_f(&mpz, &mpf);

        // un-negate, positively
        var neg = gmp.mpz_cmp_si(&mpz, 0) < 0;
        if (neg) {
            const words = gmp.mpz_size(&mpz);
            var limbPtr = gmp.mpz_limbs_write(&mpz, @intCast(c_long, words));
            var i : usize = 0;
            while (i < words) : (i += 1) limbPtr[i] = ~limbPtr[i];
            gmp.mpz_limbs_finish(&mpz, @intCast(gmp.mp_size_t, words));
            gmp.mpz_add_ui(&mpz, &mpz, 1);
        }

        // now get data
        var strptr = gmp.mpz_get_str(null, 16, &mpz);
        const slen = strlen(strptr);
        const istr = strptr[0..slen];
        std.debug.print("##> {}: {s}\n", .{istr.len, istr});

        var list = std.ArrayList(u8).init(alloc);
        var i : usize = 0;
        const p = slen % 16;
        if (p != 0) {
            // partial first word, pad with zeroes
            var q : isize = @mod(-@intCast(isize, p), 16);
            try list.append('.');
            while (q != 0) : (q -= 1) {
                try list.append('0');
            }
            try std.fmt.format(list.writer(), "{s}", .{istr[i..p]});
            i = p;
        }
        // full words
        while (i + 16 <= slen) : (i += 16) {
            try std.fmt.format(list.writer(), ".{s}", .{istr[i..i+16]});
        }
        if (i < slen) {
            // partial last word
            var q = slen - i;
            try std.fmt.format(list.writer(), ".{s}", .{istr[i..i + q]});
            while (q < 16) : (q += 1) {
                try list.append('0');
            }
        }

        const bistr = list.toOwnedSlice();
        // std.debug.print("ii> {s}\n", .{bistr});

        return bistr;
    }

    /// Write our hex format (.XXXX.XXX.XXX) to "normal" decimal
    /// caller owns returned array
    pub fn write_decimal(alloc: Allocator, istr: []const u8, bits: u16) ![:0]u8 {
        var mpz : gmp.mpz_t = undefined;
        gmp.mpz_init(&mpz);
        defer gmp.mpz_clear(&mpz);

        // gmp treats negative numbers funny; need to handle this way
        const data = try alloc.alloc(u64, istr.len / 17); defer alloc.free(data);
        var p : usize = 0;
        var words : c_uint = 0;
        while (p < istr.len) : (p += 17) {
            var w : u64 = try std.fmt.parseInt(u64, istr[p+1..p+17], 16);
            // gmp.mpz_mul_2exp(&mpz, &mpz, 64);
            // gmp.mpz_add_ui(&mpz, &mpz, w);
            data[words] = w;
            words += 1;
        }
        // std.debug.print("{} - {}\n", .{p, istr.len});
        if (p != istr.len) unreachable;

        var limbPtr = gmp.mpz_limbs_write(&mpz, words);
        // negative?
        if (data[0] >= 0x8000_0000_0000_0000) {
            var i : usize = 0;
            while (i < words) : (i +=1) limbPtr[words - i - 1] = ~data[i];
            gmp.mpz_limbs_finish(&mpz, -@intCast(gmp.mp_size_t, words));
            gmp.mpz_sub_ui(&mpz, &mpz, 1);
        } else {
            var i : usize = 0;
            while (i < words) : (i +=1) limbPtr[words - i - 1] = data[i];
            gmp.mpz_limbs_finish(&mpz, @intCast(gmp.mp_size_t, words));
        }

        // make a float which is the exact large integer
        var mpf : gmp.mpf_t = undefined;
        gmp.mpf_init(&mpf);
        defer gmp.mpf_clear(&mpf);
        gmp.mpf_set_prec(&mpf, (words << 6) + 8);

        // get the int into the float
        gmp.mpf_set_z(&mpf, &mpz);

        // divide down to a float
        gmp.mpf_div_2exp(&mpf, &mpf, (words << 6) - 8);

        // truncate to desired bits
        gmp.mpf_set_prec(&mpf, bits + 8);

        var exp : gmp.mp_exp_t = undefined;
        var strptr = gmp.mpf_get_str(null, &exp, 10, 256, &mpf);
        const slen = strlen(strptr);
        // std.debug.print("xx> ({}){s}\n", .{exp, strptr[0..slen]});

        var list = std.ArrayList(u8).init(alloc);
        var i : usize = 0;
        if (strptr[0] == '-') {
            try list.append('-');
            i += 1;
        }
        if (exp < 0) {
            try list.append('.');
            while (exp < 0) : (exp += 1) {
                try list.append('0');
            }
            try list.appendSlice(strptr[i..slen]);
        } else if (exp == 0) {
            try list.append('0');
            try list.append('.');
            try list.appendSlice(strptr[i..slen]);
        } else {
            try list.appendSlice(strptr[i..i+@intCast(usize, exp)]);
            try list.append('.');
            try list.appendSlice(strptr[i+@intCast(usize, exp)..slen]);
            while (exp > 0) : (exp -= 1) {
                try list.append('0');
            }
        }

        const bistr = try list.toOwnedSliceSentinel(0);
        // std.debug.print("dd> {s}\n", .{bistr});
        return bistr;
    }
};


test "file hexnums" {
    const alloc = std.testing.allocator;
    {
        //           0123456789abcdef
        const in = ".0000800000000000";
        var dec = try DecimalMath.write_decimal(alloc, in, 64); defer alloc.free(dec);
        try testing.expectEqualStrings(".001953125", dec);
        var hex = try DecimalMath.parse_decimal(alloc, dec, 64); defer alloc.free(hex);
        try testing.expectEqualStrings(in, hex);
    }
    {
        //           0123456789abcdef 0123456789abcdef
        const in = ".0000800000000000.0000000000000000";
        // more bits than expected
        var dec = try DecimalMath.write_decimal(alloc, in, 64); defer alloc.free(dec);
        try testing.expectEqualStrings(".001953125", dec);
        var hex = try DecimalMath.parse_decimal(alloc, dec, 128); defer alloc.free(hex);
        try testing.expectEqualStrings(in, hex);
    }
    {
        //           0123456789abcdef 0123456789abcdef
        const in = ".0000800000000000.0000000000000008";
        var dec = try DecimalMath.write_decimal(alloc, in, 64); defer alloc.free(dec);
        try testing.expectEqualStrings(".001953125000000000000000000000000006018531", dec);
        // clip precision
        var hex = try DecimalMath.parse_decimal(alloc, dec, 64); defer alloc.free(hex);
        try testing.expectEqualStrings(".0000800000000000", hex);
    }
    {
        //           0123456789abcdef
        const in = ".0100000000000000";
        var dec = try DecimalMath.write_decimal(alloc, in, 64); defer alloc.free(dec);
        try testing.expectEqualStrings("1.0", dec);
        var hex = try DecimalMath.parse_decimal(alloc, dec, 64); defer alloc.free(hex);
        try testing.expectEqualStrings(".0100000000000000", hex);
    }
    {
        //           0123456789abcdef
        const in = ".ff00000000000000";
        var dec = try DecimalMath.write_decimal(alloc, in, 64); defer alloc.free(dec);
        try testing.expectEqualStrings("-1.0", dec);
        var hex = try DecimalMath.parse_decimal(alloc, dec, 64); defer alloc.free(hex);
        try testing.expectEqualStrings(".ff00000000000000", hex);
    }
}

test "file decnums" {
    const alloc = std.testing.allocator;
    {
        const v : f64 = (1.0 / 1048576.0) * 123456;
        var in = try std.fmt.allocPrintZ(alloc, "{d:16}", .{v}); defer alloc.free(in);
        std.debug.print("--> {s}\n", .{in});
        var hex = try DecimalMath.parse_decimal(alloc, in, 64); defer alloc.free(hex);
        try testing.expectEqualStrings(".001e240000000000", hex);
        var dec = try DecimalMath.write_decimal(alloc, hex, 64); defer alloc.free(dec);
        try testing.expectEqualStrings(in, dec);
    }
    {
        const v : f64 = (1.0 / 1048576.0 / 65536.0) * 123456789012;
        // 1.7965327280689962
        var in = try std.fmt.allocPrintZ(alloc, "{d:16}", .{v}); defer alloc.free(in);
        std.debug.print("--> {s}\n", .{in});
        var hex = try DecimalMath.parse_decimal(alloc, in, 64); defer alloc.free(hex);
        try testing.expectEqualStrings(".01cbe991a13ffffe", hex);
        var dec = try DecimalMath.write_decimal(alloc, hex, 64); defer alloc.free(dec);
        // too precise, so extra digits
        try testing.expectEqualStrings("1.7965327280689961930715270455038989894090", dec);

        var dec2 = try DecimalMath.write_decimal(alloc, hex, 128); defer alloc.free(dec2);
        // a little better
        try testing.expectEqualStrings("1.79653272806899619307152704550389898940920829772949218750", dec2);

        // but this gets real precision
        var hex2 = try DecimalMath.parse_decimal(alloc, in, 128); defer alloc.free(hex2);
        try testing.expectEqualStrings(".01cbe991a13ffffe.7fcec9d85e809502", hex2);
        var dec3 = try DecimalMath.write_decimal(alloc, hex2, 128); defer alloc.free(dec3);
        // best we can get
        try testing.expectEqualStrings("1.7965327280689961999999999999999999995532760065860857476430", dec3);

    }
    {
        // 174   |                                                                                                                                                                              |
        var in = "0.39001807130093577831202936299725061179962432297641336215752898092430021221605182748563970878207227694639197294204891273579776187087084006441994806701999929579553600913981309347876057028568242043043158561091208702165728233337881810870731858475214098927501068403786243207505510898353565649399508688234219336575480457759195484856559999342027284671912897720236968725844641511183359441099508188188002174200656217974760066339717742919140401845770741563017036784716375173286823097029696327808778733015060424804687500000";
        var hex = try DecimalMath.parse_decimal(alloc, in, 576); defer alloc.free(hex);
        try testing.expectEqualStrings(".0063d8396d1625de.94b11d3317e7226c.a86153fc8b75e82e.964623c31ce9c8a8.aa74f16ede5c476c.14e747ecf1b2ae87.ec3905dae575b822.72aab22ff6050e64.4af0da41acbbae03", hex);
        var dec = try DecimalMath.write_decimal(alloc, hex, 576); defer alloc.free(dec);
        try testing.expectEqualStrings("0.39001807130093577831202936299725061179962432297641336215752898092430021221605182748563970878207227694639197294204891273579776187087084006441994806701999929579553600913981", dec[0..172]);

    }
}

///////////

pub const RenderedFile = struct {
    const Self = @This();

    params: *MandelParams,
    storage : *MandelStorage,
    alloc: Allocator,

    pub fn init(params: *MandelParams, storage: *MandelStorage) !Self {
        return Self{
            .params = params,
            .storage = storage,
            // .alloc = arena_instance.allocator()
            .alloc = std.heap.page_allocator,
        };
    }

    pub fn deinit(self: *Self) void {
        _ = self;
    }

    pub fn load(self: *Self, path: []const u8) !void {
        var file = try std.fs.cwd().openFile(path, .{});
        defer file.close();

        var reader = file.reader();
        var version = try self.read_string(reader); defer self.alloc.free(version);
        if (!std.mem.eql(u8, version, "jmandel_1.1")) {
            return error.UnsupportedVersion;
        }

        try self.load_version_1_1(reader, (try file.stat()).size);
    }

    fn load_version_1_1(self: *Self, reader: File.Reader, fileSize: u64) !void {
        var cxStr = try self.read_string(reader); defer self.alloc.free(cxStr);
        var cyStr = try self.read_string(reader); defer self.alloc.free(cyStr);

        {
            const zoomStr = try self.read_string(reader); defer self.alloc.free(zoomStr);
            self.params.zoom = try std.fmt.parseInt(u16, zoomStr, 10);
        }

        var bits = self.params.getDefaultIntSize() << 6;
        self.params.words = bits >> 6;
        self.params.setCx(self.storage.alloc, try DecimalMath.parse_decimal(self.storage.alloc, cxStr, config.MAXBITS));
        self.params.setCy(self.storage.alloc, try DecimalMath.parse_decimal(self.storage.alloc, cyStr, config.MAXBITS));

        {
            const minItersStr = try self.read_string(reader); defer self.alloc.free(minItersStr);
            var minIters = try std.fmt.parseInt(u32, minItersStr, 10);
            _ = minIters;
        }
        {
            const maxItersStr = try self.read_string(reader); defer self.alloc.free(maxItersStr);
            var maxIters = try std.fmt.parseInt(u32, maxItersStr, 10);
            _ = maxIters;
        }

        const psxStr = try self.read_string(reader); defer self.alloc.free(psxStr);
        const psyStr = try self.read_string(reader); defer self.alloc.free(psyStr);
        const psx = try std.fmt.parseInt(u32, psxStr, 10);
        const psy = try std.fmt.parseInt(u32, psyStr, 10);
        self.params.sx = psx;
        self.params.sy = psy;

        const exp : u32 = psx * psy;
        var iters = try self.alloc.alloc(i32, exp); defer self.alloc.free(iters);
        std.mem.set(i32, iters, 0);
        try self.decompress_iters(iters, reader, fileSize);

        try switch (self.params.words) {
            inline 1...bignum.MAXWORDS => |w| self.decompress_blocks(bignum.BigInt(w << 6), psx, psy, iters),
            else => unreachable,
        };

        std.debug.print("Loaded blocks\n", .{});
    }

    fn decompress_iters(self: Self, iters : []i32, reader : File.Reader, fileSize: u64) !void {
        const zmem = try self.alloc.alloc(u8, fileSize); defer self.alloc.free(zmem);
        const zlen = try reader.readAll(zmem);

        var zstr : zlib.z_stream = undefined;
        std.mem.set(u8, @ptrCast([*]u8, &zstr)[0..@sizeOf(zlib.z_stream)], 0);
        // window size (15, default) gzip header (+32)
        if (zlib.inflateInit2(&zstr, 15 + 32) != zlib.Z_OK) return error.ZLibError;
        zstr.next_in = zmem.ptr;
        zstr.avail_in = @intCast(c_uint, zlen);

        var i8ptr : usize = 0;
        var iters8 = @ptrCast([*]u8, iters);

        while (zstr.avail_in > 0) {
            std.debug.print("zlib iter {}...\n", .{ zstr.avail_in });
            const avail = @intCast(c_uint, iters.len * 4 - i8ptr);
            zstr.avail_out = avail;
            zstr.next_out = iters8 + i8ptr;

            const res = zlib.inflate(&zstr, if (zstr.avail_in == 0) zlib.Z_FINISH else zlib.Z_NO_FLUSH);
            if (res == zlib.Z_STREAM_END or zstr.avail_out == 0) break;
            if (res == zlib.Z_DATA_ERROR) {
                std.debug.print("zlib error: {s}\n", .{zstr.msg});
                return; // error.FileFormat;
            }

            i8ptr = zstr.total_out;
        }

        _ = zlib.inflateEnd(&zstr);
    }

    fn decompress_blocks(self : *Self, comptime BigIntType : type, psx : u32, psy : u32, iters: []const i32) !void {
        _ = psy;
        var workLoad = try MandelParams.BlockWorkMaker(BigIntType).init(self.storage.alloc, self.params.*, self.storage, self.params.sx, 0, 0);
        defer workLoad.deinit(self.storage.alloc);

        std.debug.print("iters={}, psx={}, psy={}\n", .{ iters.len, self.params.sx, self.params.sy });

        const ps = self.storage.blockSize;
        var c : u32 = 0;
        for (workLoad.blocks) |block| {
            const px = block.px;
            const py = block.py;

            var oy : u32 = 0;
            while (oy < ps) : (oy += 1) {
                var ox : u32 = 0;
                while (ox < ps) : (ox += 1) {
                    const iterBE : i32 = iters[(py + oy) * psx + (px + ox)];
                    const iter = @byteSwap(iterBE);
                    block.data.setIter(ox, oy, iter);
                    c += 1;

                    // if (iter != l) {
                    //     std.debug.print("{} ", .{iter});
                    //     l = iter;
                    //     i += 1;
                    //     if (@rem(i, 16) == 0) std.debug.print("\n", .{});
                    // }

                }
            }
        }
    }

    fn read_string(self: Self, reader : File.Reader) ![:0]u8 {
        var length = try reader.readIntBig(u16);
        var str = try self.alloc.allocSentinel(u8, length, 0);
        var len = try reader.read(str);
        if (len != length) return error.BadFileFormat;
        std.debug.print("--> {s}\n", .{str});
        return str;
    }

    pub fn save(self: *Self, path: []const u8) !void {
        var file = try std.fs.cwd().createFile(path, .{.truncate = true});
        defer file.close();

        var writer = file.writer();
        try self.write_string_free(writer, try self.alloc.dupe(u8, "jmandel_1.1"));

        const iters = try self.alloc.alloc(i32, self.params.sx * self.params.sy); defer self.alloc.free(iters);
        std.mem.set(i32, iters, 0);

        var minIters : u32 = undefined;
        var maxIters : u32 = undefined;

        try switch (self.params.words) {
            inline 1...bignum.MAXWORDS => |w| self.compress_blocks(bignum.BigInt(w << 6), self.params.sx, self.params.sy, iters, &minIters, &maxIters),
            else => unreachable,
        };

        const cxStr = try DecimalMath.write_decimal(self.alloc, self.params.cx, config.MAXBITS);
        const cyStr = try DecimalMath.write_decimal(self.alloc, self.params.cy, config.MAXBITS);
        try self.write_string_free(writer, cxStr);
        try self.write_string_free(writer, cyStr);

        try self.write_string_free(writer, try std.fmt.allocPrint(self.alloc, "{}", .{ self.params.zoom }));

        try self.write_string_free(writer, try std.fmt.allocPrint(self.alloc, "{}", .{ minIters }));
        try self.write_string_free(writer, try std.fmt.allocPrint(self.alloc, "{}", .{ maxIters }));

        try self.write_string_free(writer, try std.fmt.allocPrint(self.alloc, "{}", .{ self.params.sx }));
        try self.write_string_free(writer, try std.fmt.allocPrint(self.alloc, "{}", .{ self.params.sy }));

        try self.compress_iters(iters, writer);
    }

    fn write_string_free(self: Self, writer : File.Writer, str: []u8) !void {
        try writer.writeIntBig(u16, @intCast(u16, str.len));
        try writer.writeAll(str);
        self.alloc.free(str);
    }

    fn compress_blocks(self : *Self, comptime BigIntType : type, psx : u32, psy : u32, iters: []i32, minIters: *u32, maxIters: *u32) !void {
        _ = psy;
        var workLoad = try MandelParams.BlockWorkMaker(BigIntType).init(self.storage.alloc, self.params.*, self.storage, self.params.sx, 0, 0);
        defer workLoad.deinit(self.storage.alloc);

        std.debug.print("iters={}, psx={}, psy={}\n", .{ iters.len, self.params.sx, self.params.sy });

        minIters.* = std.math.maxInt(u32);
        maxIters.* = 0;

        const ps = self.storage.blockSize;
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

                    if (iter > 0) {
                        minIters.* = std.math.min(minIters.*, iter);
                        maxIters.* = std.math.max(maxIters.*, @intCast(u32, iter));
                    }

                    c += 1;
                }
            }
        }
    }

    fn compress_iters(self: Self, iters: []i32, writer: File.Writer) !void {
        var obuf = try self.alloc.alloc(u8, 4096); defer self.alloc.free(obuf);

        var zstr : zlib.z_stream = undefined;
        std.mem.set(u8, @ptrCast([*]u8, &zstr)[0..@sizeOf(zlib.z_stream)], 0);

        var i8len : usize = iters.len * 4;

        zstr.next_in = @ptrCast([*] u8, iters);
        zstr.avail_in = @intCast(c_uint, i8len);

        const ires = zlib.deflateInit2(&zstr,
            zlib.Z_DEFAULT_COMPRESSION, zlib.Z_DEFLATED,
            15 + 16, // window size (15, default) gzip header (+16)
            zlib.MAX_MEM_LEVEL, // memLevel
            zlib.Z_DEFAULT_STRATEGY);
        if (ires != zlib.Z_OK) {
            std.debug.print("init error: {}\n", .{ires});
            return error.ZLibError;
        }

        zstr.next_in = @ptrCast([*] u8, iters);
        zstr.avail_in = @intCast(c_uint, i8len);

        while (true) {
            zstr.avail_out = @intCast(c_uint, obuf.len);
            zstr.next_out = obuf.ptr;

            // const res = zlib.deflate(&zstr, if (zstr.avail_in == 0) zlib.Z_FINISH else zlib.Z_NO_FLUSH);
            const res = zlib.deflate(&zstr, zlib.Z_FINISH);
            if (res != zlib.Z_STREAM_END and res != zlib.Z_OK) {
                if (zstr.msg != null)
                    std.debug.print("zlib error: {s}\n", .{zstr.msg})
                else
                    std.debug.print("zlib error: {}\n", .{res});
                return error.CompressError;
            }

            const olen = obuf.len - zstr.avail_out;
            // std.debug.print("olen={}\n", .{olen});
            _ = try writer.writeAll(obuf[0..olen]);

            if (zstr.avail_out != 0) break;
        }

        _ = zlib.deflateEnd(&zstr);
    }
};

/////////////////

fn run_test(comptime callable: anytype) !void {
    const alloc = std.testing.allocator;
    // const alloc = std.heap.page_allocator;
    var storage = try MandelStorage.init(alloc, 128);
    defer storage.deinit();
    var params = try MandelParams.init(alloc);
    defer params.deinit(alloc);

    try callable.doit(&storage, &params);

}


test "file load1" {
    try run_test(struct{
        pub fn doit(storage: *MandelStorage, params: *MandelParams) !void {
            var file = try RenderedFile.init(params, storage);
            defer file.deinit();
            try file.load("data/zoom10_000.dat");
        }
    });
}
test "file save1" {
    try run_test(struct{
        pub fn doit(storage: *MandelStorage, params: *MandelParams) !void {
            var file = try RenderedFile.init(params, storage);
            defer file.deinit();
            try file.load("data/zoom10_400.dat");
            try file.save("data/zoom10_400.dat.bak");
        }
    });
}
