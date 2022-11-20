const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;
const json = std.json;

const config = @import("config.zig");
const BigIntMax = config.BigIntMax;

const bignum = @import("bignum.zig");
const BigInt = bignum.BigInt;
const BigFixedFloat = bignum.BigFixedFloat;

const algo = @import("algo.zig");
const MandelAlgo = algo.Mandel;

const XY = @import("types.zig").XY;

pub const BigCoord = struct {
    const Self = @This();

    x: BigIntMax,
    y: BigIntMax,

    pub fn init(xi : BigIntMax, yi: BigIntMax) Self {
        return .{ .x = xi, .y = yi };
    }
    pub fn initFrom(comptime T : type, x : T, y: T) Self {
        const shiftBits = config.MAXBITS - @typeInfo(T).Int.bits;
        return .{ .x = @intCast(BigIntMax, x) << shiftBits, .y = @intCast(BigIntMax, y) << shiftBits };
    }

    pub fn hash(self: Self) u64 {
        var hasher = std.hash.Wyhash.init(0);
        std.hash.autoHashStrat(&hasher, self.x, .Deep);
        std.hash.autoHashStrat(&hasher, self.y, .Deep);
        const h = hasher.final();
        return h;
    }
    pub fn eql(self: BigCoord, other: BigCoord) bool {
        return self.x == other.x and self.y == other.y;
    }

    pub fn to_string(self: *const BigCoord, alloc: Allocator) ![]u8 {
        var xs = try bignum.printBig(alloc, BigIntMax, self.x); defer alloc.free(xs);
        var ys = try bignum.printBig(alloc, BigIntMax, self.y); defer alloc.free(ys);
        return try std.fmt.allocPrint(alloc, "[{s}, {s}]", .{ xs, ys });
    }
};

/// Calculated block data, meaningful relative to some Params
pub const BlockData = struct {
    const Self = @This();

    /// size of box in pixels
    sz: u16,
    /// size of box in shift
    szsh: u5,

    /// iter counts (0 if not calculated, +ve = iter count, -ve = last failed count)
    iters: []i32,

    pub fn init(self: *Self, alloc: Allocator, sz: u16) !void {
        self.sz = sz;
        self.szsh = std.math.log2_int(u16, sz);
        self.iters = try alloc.alloc(i32, @intCast(u32, sz) * sz);
        std.mem.set(i32, self.iters, 0);
    }

    pub fn deinit(self: *Self, arena: Allocator) void {
        arena.free(self.iters);
    }

    pub inline fn iter(self: *const Self, x: u32, y: u32) i32 {
        return self.iters[(y << self.szsh) + x];
    }

    pub inline fn setIter(self: *Self, x: u32, y: u32, it: i32) void {
        self.iters[(y << self.szsh) + x] = it;
    }
};

/// All the data calculated at a given zoom level
pub const MandelLayer = struct {
    const Self = @This();

    // we hold pointers to BlockData so the memory in the hashmap doesn't
    // move its contents for us...
    const BlockMap = std.HashMap(BigCoord, *BlockData, struct {
        pub fn hash(_: anytype, a: BigCoord) u64 { return a.hash(); }
        pub fn eql(_: anytype, a: BigCoord, b: BigCoord) bool { return a.eql(b); }
    }, 80);

    arena: Allocator,
    zoom: u16,
    /// block size (pixels)
    blockSize: u16,
    /// cache of blocks
    blocks: BlockMap,

    pub fn init(self: *Self, zoom: u16, blockSize: u16) !void {
        // self.arena = arena_instance.allocator();
        self.arena = std.heap.page_allocator;
        self.zoom = zoom;
        self.blockSize = blockSize;
        self.blocks = BlockMap.init(self.arena);
    }

    pub fn deinit(self: *Self) void {
        while (true) {
            var it = self.blocks.iterator();
            var ent = it.next();
            if (ent == null) break;
            ent.?.value_ptr.*.deinit(self.arena);
            _ = self.blocks.remove(ent.?.key_ptr.*);
        }
        self.blocks.clearAndFree();
        self.blocks.deinit();
    }

    /// Get the block
    pub fn get(self: *const Self, coord: BigCoord) ?*BlockData {
        return self.blocks.get(coord);
    }

    /// Get the block
    pub fn ensure(self: *Self, coord: BigCoord) !*BlockData {
        {
            var blockPtr = self.get(coord);
            if (blockPtr != null) {
                return blockPtr.?;
            }
        }
        {
            var cs = try coord.to_string(self.arena); defer self.arena.free(cs);
            var res = try self.blocks.getOrPut(coord);
            var blockPtr = res.value_ptr;
            if (!res.found_existing) {
                blockPtr.* = try self.arena.create(BlockData);
                try blockPtr.*.init(self.arena, self.blockSize);
            }
            return blockPtr.*;
        }
    }
};

/// All the data we've computed
pub const MandelStorage = struct {
    const Self = @This();

    // we hold pointers to MandelLayer so the memory in the hashmap doesn't
    // move its contents for us (or copy it...)...
    const LayerMap = std.HashMap(u16, *MandelLayer, struct {
        pub fn hash(_: anytype, a: u16) u64 { return a; }
        pub fn eql(_: anytype, a: u16, b: u16) bool { return a == b; }
    }, 80);

    alloc: Allocator,
    blockSize: u16,
    layers: LayerMap,

    pub fn init(alloc: Allocator, blockSize: u16) !Self {
        return .{
            .alloc = alloc,
            .blockSize = blockSize,
            .layers = LayerMap.init(alloc)
        };
    }

    pub fn deinit(self: *Self) void {
        // var it = self.layers.iterator();
        // while (it.next()) |ent| {
        //     ent.value_ptr.*.deinit();
        //     // self.alloc.destroy(ent.value_ptr);
        // }
        while (true) {
            var it = self.layers.iterator();
            var ent = it.next();
            if (ent == null) break;
            ent.?.value_ptr.*.deinit();
            self.alloc.destroy(ent.?.value_ptr.*);
            _ = self.layers.remove(ent.?.key_ptr.*);
        }

        self.layers.clearAndFree();
        self.layers.deinit();
    }

    pub fn get(self: *Self, zoom: u16) ?*const MandelLayer {
        return self.layers.get(zoom);
    }

    pub fn ensure(self: *Self, zoom: u16, sz: u16) !*MandelLayer {
        var gop = try self.layers.getOrPut(zoom);
        var layer = gop.value_ptr;
        if (!gop.found_existing) {
            std.debug.print("new zoom level layer {}\n", .{zoom});
            layer.* = try self.alloc.create(MandelLayer);
            try layer.*.init(zoom, sz);
        }
        return layer.*;
    }

    pub fn remove(self: *Self, zoom: u16) void {
        var layer = self.layers.get(zoom);
        if (layer != null) {
            _ = self.layers.remove(zoom);
            layer.?.deinit();
            self.alloc.destroy(layer.?);
        }
    }

};


test "bigcoord" {
    const c : BigCoord = BigCoord.init(0x123456, 0xabcdef);
    const alloc = std.testing.allocator;

    var cs = try c.to_string(alloc); defer alloc.free(cs);
    std.debug.print("{s}\n", .{cs});
}
