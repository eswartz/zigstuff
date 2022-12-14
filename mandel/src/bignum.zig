const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;

const config = @import("config.zig");

pub const MAXWORDS = config.MAXWORDS;

pub fn BigInt(comptime size: u16) type {
    return @Type(.{ .Int = .{ .signedness = .signed, .bits = size } });
}

pub const NativeFloat = f128;
pub const NativeInt = i64;


fn faddShift_T(comptime T: type, alloc: Allocator, cstr: *[]u8, addend: NativeFloat, shift: u6) !void {
    const f = try parseBig(T, cstr.*);
    const fadd: T = BigFixedFloat(T, 8).fromFloat(addend);
    // var fstr = try printBig(alloc, T, fadd);
    // defer alloc.free(fstr);
    // std.debug.print("from {s} fadd={s}, shift={}\n", .{ cstr, fstr, shift });
    const faddend: T = fadd >> shift;
    // var cstr2 = try printBig(alloc, T, faddend);
    // defer alloc.free(cstr2);

    const res = try printBig(alloc, T, f + faddend);
    // std.debug.print("+ {s} = res={s}\n", .{ cstr2, res });
    alloc.free(cstr.*);
    cstr.* = res;
}

pub fn faddShift(alloc: Allocator, cstr: *[]u8, words: u16, addend: NativeFloat, shift: u6) !void {
    return switch (words) {
        inline 1...MAXWORDS => |v| faddShift_T(BigInt(64 * v), alloc, cstr, addend, shift),
        else => unreachable,
    };
}

fn fadj_T(comptime BigIntType: type, alloc: Allocator, cstr: *[]u8, span: NativeFloat, offs: i32, range: i32) !void {
    const f = try parseBig(BigIntType, cstr.*);

    const fspan: BigIntType = BigFixedFloat(BigIntType, 8).fromFloat(span);

    const shift = std.math.log2_int(u32, @intCast(u32, range));

    // var fstr = try printBig(alloc, BigIntType, fspan);
    // defer alloc.free(fstr);

    // std.debug.print("from {s} span={s}, offs={}, range={}, shift={}\n", .{ cstr, fstr, offs, range, shift });

    // params.*.cx += @intToFloat(NativeFloat, me.x - @intCast(i32, params.*.sx >> 1)) * params.*.span() / @intToFloat(NativeFloat, params.*.sx);
    const BigFloat2 = BigInt(@typeInfo(BigIntType).Int.bits * 2);
    const pos2 = (@as(BigFloat2, offs) * @as(BigFloat2, fspan)) >> shift;
    const pos = @intCast(BigIntType, pos2) + f;

    // var cstr2 = try printBig(alloc, BigIntType, pos);
    // defer alloc.free(cstr2);

    const res = try printBig(alloc, BigIntType, pos);
    // std.debug.print("+ {s} = res={s}\n", .{ cstr2, res });

    alloc.free(cstr.*);
    cstr.* = res;
}

pub fn fadj(alloc: Allocator, cstr: *[]u8, words: u16, span: NativeFloat, offs: i32, range: i32) !void {
    return switch (words) {
        inline 1...MAXWORDS => |v| fadj_T(BigInt(64 * v), alloc, cstr, span, (offs) - (range >> 1), range),
        else => unreachable,
    };
}

fn falign_T(comptime BigIntType: type, alloc: Allocator, cstr: *[]u8) !void {
    const f = try parseBig(BigIntType, cstr.*);
    const res = try printBig(alloc, BigIntType, f);
    alloc.free(cstr.*);
    cstr.* = res;
}

pub fn falign(alloc: Allocator, cstr: *[]u8, words: u16) !void {
    return switch (words) {
        inline 1...MAXWORDS => |v| falign_T(BigInt(64 * v), alloc, cstr),
        else => unreachable,
    };
}


fn hexval(c: u8) u32 {
    return if (c >= 'a') c - 'a' + 10 else if (c >= 'A') c - 'A' + 10 else c - '0';
}

fn parseHexFloat(comptime T: type, s: []const u8) !T {
    if (s.len == 0) {
        return error.InvalidCharacter;
    }

    var i: usize = 0;
    const negative = s[i] == '-';
    if (s[i] == '-' or s[i] == '+') {
        i += 1;
    }
    if (s.len == i) {
        return error.InvalidCharacter;
    }

    // skip '0x'
    if (s[i] != '0') unreachable;
    i += 1;
    if (s[i] != 'x') unreachable;
    i += 1;

    // int part -- always one
    var f: T = 0;

    if (s[i] != '1') unreachable;
    i += 1;
    if (s[i] != '.') unreachable;
    i += 1;

    // mantissa
    var ep: i32 = 0;
    var fi: i128 = 1;
    while (i < s.len) : (i += 1) {
        if (s[i] == 'p') {
            i += 1;
            break;
        }
        fi = fi * 16 + hexval(s[i]);
        ep += 4;
        // std.debug.print("fi = {}\n", .{ fi });
    }
    f = std.math.scalbn(@intToFloat(NativeFloat, fi), -ep);
    // std.debug.print("f = {}, ep = {}\n", .{ f, ep });

    // exponent
    var ex: i32 = 0;
    if (s[i] != '-') unreachable;
    i += 1;

    // std.debug.print("exp = {s}\n", .{ s[i..] });
    while (i < s.len) : (i += 1) {
        ex = ex * 10 + (s[i] - '0');
    }
    // std.debug.print("ex = {}, ep = {}\n", .{ ex, ep });

    f *= std.math.scalbn(@floatCast(NativeFloat, 1.0), -ex);

    return if (negative) -f else f;
}

pub fn parseBig(comptime IntType: type, str: []const u8) !IntType {
    var v: IntType = 0;
    const wordCount = @typeInfo(IntType).Int.bits >> 6;
    var words = @ptrCast([*]u64, &v)[0..wordCount];

    var idx: u32 = 0;
    var vidx: u32 = 0;
    while (idx < str.len and vidx < wordCount) {
        if (str[idx] != '.') return error.InvalidCharacter;
        idx += 1;
        words[wordCount - vidx - 1] = try std.fmt.parseInt(u64, str[idx .. idx + 16], 16);
        idx += 16;
        vidx += 1;
    }

    return v;
}

pub fn printBigArray(alloc: Allocator, words: []const u64) ![]u8 {
    var list = std.ArrayList(u8).init(alloc);

    const wordCount = words.len;

    var vidx: u32 = 0;
    while (vidx < wordCount) : (vidx += 1) {
        try std.fmt.format(list.writer(), ".{x:0>16}", .{words[wordCount - vidx - 1]});
    }
    return list.toOwnedSlice();
}

/// Returns reference to slice, not copied (hence, must be inlined)
pub inline fn intWords(comptime IntType: type, x: IntType) []const u64 {
    const wordCount = @typeInfo(IntType).Int.bits >> 6;
    const words = @ptrCast([*]const u64, &x)[0..wordCount];
    return words;
}

fn intsAdd(comptime IntType: type, a: []const u64, b: []const u64, s: []u64) void {
    var ai = @ptrCast(*const IntType, a.ptr);
    var bi = @ptrCast(*const IntType, b.ptr);
    @ptrCast(*IntType, s.ptr).* = ai.* + bi.*;
}

pub fn intWordsAdd(alloc: Allocator, x: []const u64, y: []const u64) ![]u64 {
    var s = try alloc.alloc(u64, x.len);
    switch(x.len) {
        inline 1...16 => |v| intsAdd(BigInt(v * 64), x, y, s),
        else => unreachable,
    }
    return s;
}


pub fn printBig(alloc: Allocator, comptime IntType: type, x: IntType) ![]u8 {
    return try printBigArray(alloc, intWords(IntType, x));
}

pub fn BigFixedFloat(comptime BigFloatType: type, comptime intBits: u16) type {
    // _ = mantissaBits;
    return struct {
        pub const Repr = BigFloatType;
        pub const mantissaBits = @typeInfo(BigFloatType).Int.bits - intBits;
        /// Convert 'v' to a big float.
        /// Unfortunately some bigint/float operations aren't handled yet,
        /// so we need to do this manually
        pub fn fromFloat(v: NativeFloat) Repr {
            // whoops, LLVM ERROR...
            // return @floatToInt(BigFloatType, av);

            // const fval = v * @intToFloat(F, (@as(BigFloatType, 1) << mantissaBits));
            // return @floatToInt(BigFloatType, fval);

            // var base : F = std.math.scalbn(@as(F, 1.0), mantissaBits);
            // var rough : i64 = @floatToInt(i64, v * base);
            // std.debug.print("rough = {x}\n", .{ rough });
            // var f : BigFloatType = @as(BigFloatType, rough) << (@typeInfo(BigFloatType).Int.bits - 64);
            // return f;

            var av: NativeFloat = std.math.fabs(v);
            const fexp = std.math.frexp(av);
            _ = fexp;
            var f: Repr = @as(Repr, @floatToInt(NativeInt, av));
            var i: isize = mantissaBits;
            while (i > 0) : (i -= 1) {
                av -= @trunc(av);
                av *= 2.0;
                f = if (av >= 1.0) (f << 1) | 1 else f << 1;
            }
            return if (v < 0) -f else f;
        }
    };
}

// pub var p : BigInt256 = 1;

// pub fn main() !void {
//     p *= 9;
//     _ = p;
// }

fn testBigEqu(comptime T: type, x: T, y: T) !void {
    const print = @import("std").debug.print;

    const alloc = std.testing.allocator;
    var s1 = try printBig(alloc, T, x); defer alloc.free(s1);
    var s2 = try printBig(alloc, T, y); defer alloc.free(s2);

    if (x == y) return;

    print("expected:\n\t{s},\ngot:\n\t{s}\n", .{ s1, s2 });
    try testing.expect(false);
}

test "basic64" {
    const BigInt64 = BigInt(64);
    const m = BigFixedFloat(BigInt64, 8);
    var x = m.fromFloat(0.5);
    var y = m.fromFloat(-0.25);
    try testBigEqu(BigInt64, @as(BigInt64, 0x0080_0000_0000_0000), x);
    try testBigEqu(BigInt64, @as(BigInt64, -0x0040_0000_0000_0000), y);
}

test "basic256" {
    const BigInt256 = BigInt(256);
    const m = BigFixedFloat(BigInt256, 8);
    var x = m.fromFloat(0.5);
    var y = m.fromFloat(-0.25);
    try testBigEqu(BigInt256, @as(BigInt256, 0x0080) << (256 - 16), x);
    try testBigEqu(BigInt256, @as(BigInt256, -0x0040) << (256 - 16), y);
    try testBigEqu(BigInt256, @as(
        BigInt256,
        -0x00555555_55555555_55555555_55555540_00000000_00000000_00000000_00000000,
    ), m.fromFloat(-1.0 / 3.0));
}

test "repr256" {
    const BigInt256 = BigInt(256);
    const m = BigFixedFloat(BigInt256, 8);
    var f: NativeFloat = 0.5;
    var exp = m.fromFloat(f);
    var i: usize = 0;
    while (i < 256 - 16) : (i += 1) {
        var x = m.fromFloat(f);
        try testBigEqu(BigInt256, exp, x);
        f /= 2.0;
        exp >>= 1;
    }
}
