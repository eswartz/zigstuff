const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;
const json = std.json;

const config = @import("config.zig");
const BigIntMax = config.BigIntMax;

const bignum = @import("bignum.zig");
const Mandel = bignum.Mandel;

const boxfill = @import("boxfill.zig");
const XY = boxfill.XY;

const rb = @import("rb.zig");

var arena_instance = std.heap.ArenaAllocator.init(std.heap.raw_c_allocator);

pub const Params = struct {
    sx: u32,
    sy: u32,
    zoom: u16,
    words: u16,
    magShift: u5,
    iters: u32,
    cx: []u8,
    cy: []u8,

    pub fn init(alloc: Allocator) !Params {
        return Params{
            .sx = 1024,
            .sy = 1024,
            .cx = try alloc.dupe(u8, ".0000000000000000"),
            .cy = try alloc.dupe(u8, ".0000000000000000"),
            .zoom = 0,
            .magShift = 2,
            .words = 1,
            .iters = 300,
        };
    }

    pub fn deinit(self: *Params, alloc: Allocator) void {
        alloc.free(self.cx);
        alloc.free(self.cy);
    }

    pub fn setCx(self: *Params, alloc: Allocator, str: []u8) void {
        alloc.free(self.cx);
        self.cx = str;
    }
    pub fn setCy(self: *Params, alloc: Allocator, str: []u8) void {
        alloc.free(self.cy);
        self.cy = str;
    }

    pub fn span(self: Params, calcSize : u32) bignum.NativeFloat {
        const fullSpan = std.math.scalbn(@floatCast(bignum.NativeFloat, 4.0), -@intCast(i16, self.zoom));
        return fullSpan / @intToFloat(bignum.NativeFloat, @divExact(self.sx, calcSize));
    }

    pub fn segments(_: Params, calcSize: u32, blockSize: u16) u32 {
        return @divExact(calcSize, blockSize);
    }

    /// Get number of words needed
    pub fn getDefaultIntSize(self: Params) u16 {
        const zoom = self.zoom;
        const size: u16 = (zoom + 64 + std.math.log2_int(u32, self.sx) + 63) >> 6;
        std.debug.print("zoom={}, isize={}\n", .{ zoom, size });
        return size;
    }

    pub fn writeConfig(self: Params, alloc: Allocator, path: []const u8) !void {
        _ = alloc;
        var file = try std.fs.cwd().createFile(path, .{ .truncate = true });
        std.debug.print("file = {}, path = {s}\n", .{ file, path });

        try json.stringify(.{
            .zoom = self.zoom,
            .iters = self.iters,
            .magShift = self.magShift,
            .cx = self.cx,
            .cy = self.cy,
            .words = self.words,
        }, .{}, file.writer());
    }

    pub fn loadConfig(self: *Params, alloc: Allocator, path: []const u8) !void {
        var file = try std.fs.cwd().openFile(path, .{});
        var parser = std.json.Parser.init(alloc, false);
        defer parser.deinit();
        // try json.stringify(self, .{}, file.reader());
        var input = try alloc.alloc(u8, 1024); defer alloc.free(input);
        var size = try file.reader().readAll(input);
        var text = input[0..size];
        std.debug.print("text = {s}\n", .{text});

        var tree = try parser.parse(text);
        self.setCx(alloc, try alloc.dupe(u8, tree.root.Object.get("cx").?.String));
        self.setCy(alloc, try alloc.dupe(u8, tree.root.Object.get("cy").?.String));
        self.zoom = @intCast(u16, tree.root.Object.get("zoom").?.Integer);
        self.iters = @intCast(u32, tree.root.Object.get("iters").?.Integer);
        self.magShift = @intCast(u5, tree.root.Object.get("magShift").?.Integer);
        const jsonBits = tree.root.Object.get("words");
        self.words = if (jsonBits) |j| @intCast(u16, j.Integer) else self.getDefaultIntSize();
        // self.*.sx = tree.root.Object.get("sx").?.Float;
        // self.*.sy = tree.root.Object.get("sy").?.Float;
    }

    pub fn BlockWorkMaker(comptime T: type) type {
        const BlockCalc = struct {
            const Self = @This();
            const BigFloat = bignum.BigFixedFloat(T, 8);

            // current data
            data: *BlockData,
            // if set, data from level zoom-1
            datam1: ?*const BlockData,
            // if set, data from level zoom+1
            datap1_00: ?*const BlockData,
            datap1_01: ?*const BlockData,
            datap1_10: ?*const BlockData,
            datap1_11: ?*const BlockData,

            // coord of box in viewport
            px: u32,
            py: u32,
            // how much to calculate
            maxIters: u32,
            // skip size from magnification
            rectSize: u32,
            // reference copy
            fxs: []T,
            // reference copy
            fys: []T,

            pub fn calculate(self: Self, stop: *std.atomic.Atomic(bool)) !bool {
                const rectSize = self.rectSize;
                const data = self.data;

                const blockSize = data.sz;
                var all : u32 = 0;
                var recalc : u32 = 0;
                var oy: u32 = 0;
                const halfBlockSize = blockSize >> 1;

                while (!stop.*.load(.Unordered) and oy < blockSize) : (oy += rectSize) {
                    const py = self.py + oy;
                    const fy = self.fys[py];
                    var ox: u32 = 0;
                    while (ox < blockSize) : (ox += rectSize) {
                        const px = self.px + ox;

                        var iters = data.iter(ox, oy);
                        if (iters == 0) {
                            if (self.datam1 != null and (ox & 1) == 0 and (oy & 1) == 0) {
                                // reuse data from outer zoom
                                iters = self.datam1.?.iter(ox / 2, oy / 2);
                            }
                            if (iters == 0) {
                                // reuse data from inner zoom?
                                var inner : ?*const BlockData =
                                    if (ox < halfBlockSize)
                                        if (oy < halfBlockSize) self.datap1_00 else self.datap1_01
                                        else if (oy < halfBlockSize) self.datap1_10 else self.datap1_11;
                                if (inner != null) {
                                    iters = inner.?.iter(((ox << 1) & (blockSize - 1)), ((oy << 1) & (blockSize - 1)));
                                }
                            }
                            data.setIter(ox, oy, iters);
                        }
                        if (iters == 0 or (iters < 0 and -iters < self.maxIters)) {
                            // need to calculate (more)
                            iters = BigFloat.calcMandelbrot(self.fxs[px], fy, self.maxIters);
                            data.setIter(ox, oy, iters);
                            recalc += 1;
                        }
                        all += 1;
                    }
                }

                if (all > 0) {
                    const perc = recalc * 16 / all;
                    std.debug.print("{s}", .{".123456789ABCDEF="[perc..perc+1]});
                }
                return recalc > 0;
            }
        };

        const Blocks = struct {
            const Self = @This();
            const BigFloatType = bignum.BigFixedFloat(T, 8);

            blocks: []BlockCalc,
            fxs: []T,
            fys: []T,

            pub fn init(alloc: Allocator, params: Params, storage: *MandelStorage,
                        calcSize : u32, xoffs : i32, yoffs : i32) !Self {
                var fxs = std.ArrayList(T).init(alloc);
                var fys = std.ArrayList(T).init(alloc);

                const pspan = params.span(calcSize);
                std.debug.print("span = {}\n", .{pspan});

                const fx0 = try bignum.parseBig(T, params.cx);
                const fy0 = try bignum.parseBig(T, params.cy);

                // split drawscape into segments
                const blockSize = storage.blockSize;
                const halfBlockSize = blockSize >> 1;
                const SEGS = params.segments(calcSize, blockSize);
                const div = @intToFloat(bignum.NativeFloat, blockSize * SEGS);
                var step = pspan / div;
                std.debug.print("segs = {}, span = {}, div = {}, step = {}\n", .{ SEGS, pspan, div, step });
                // per-pixel step
                const floatStep = BigFloatType.fromFloat(step);

                var layer = try storage.ensure(params.zoom, blockSize);
                var layerM1 = if (params.zoom > 0) storage.get(params.zoom - 1) else null;
                var layerP1 = storage.get(params.zoom + 1);

                var fss = try bignum.printBig(alloc, T, floatStep); defer alloc.free(fss);
                std.debug.print("fss = {s}\n", .{ fss });

                // calculate X and Y coords once
                {
                    const hx = @intCast(i32, calcSize / 2);
                    var px: i32 = -hx;
                    while (px < calcSize) : (px += 1) {
                        try fxs.append(fx0 + floatStep * (px + xoffs));
                    }
                }
                {
                    const hy = @intCast(i32, calcSize / 2);
                    var py: i32 = -hy;
                    while (py < calcSize) : (py += 1) {
                        try fys.append(fy0 + floatStep * (py + yoffs));
                    }
                }

                // then make blocks that reference indices into those arrays
                var blocks = std.ArrayList(BlockCalc).init(alloc);

                const rectSize: u32 = @as(u32, 1) << params.magShift;

                // get a calculation order, for quickest results
                // in the interesting area of the center
                var xys: []XY = try alloc.alloc(XY, SEGS*SEGS); defer alloc.free(xys);
                boxfill.fillBoxesInnerToOuter(xys, SEGS);

                // per block...
                for (xys) |xy| {
                    const px0 = xy.x * blockSize;
                    const py0 = xy.y * blockSize;

                    const coord = BigCoord.initFrom(T, fxs.items[px0], fys.items[py0]);
                    const blockData = try layer.ensure(coord);
                    // std.debug.print("{},{} blockData = {*}\n", .{ px0, py0, blockData });
                    const blockDataM1 = if (layerM1 != null) layerM1.?.get(coord) else null;
                    const blockDataP1_00 = if (layerP1 != null) layerP1.?.get(coord) else null;
                    const blockDataP1_10 = if (layerP1 != null) layerP1.?.get(BigCoord.initFrom(T, fxs.items[px0 + halfBlockSize], fys.items[py0])) else null;
                    const blockDataP1_01 = if (layerP1 != null) layerP1.?.get(BigCoord.initFrom(T, fxs.items[px0], fys.items[py0 + halfBlockSize])) else null;
                    const blockDataP1_11 = if (layerP1 != null) layerP1.?.get(BigCoord.initFrom(T, fxs.items[px0 + halfBlockSize], fys.items[py0 + halfBlockSize])) else null;

                    try blocks.append(.{
                        .px = px0,
                        .py = py0,
                        .maxIters = params.iters,
                        .rectSize = rectSize,
                        .fxs = fxs.items,
                        .fys = fys.items,
                        .data = blockData,
                        .datam1 = blockDataM1,
                        .datap1_00 = blockDataP1_00,
                        .datap1_01 = blockDataP1_01,
                        .datap1_10 = blockDataP1_10,
                        .datap1_11 = blockDataP1_11,
                    });
                }

                std.debug.print("Made {} blocks\n", .{blocks.items.len});
                return .{
                    .fxs = fxs.toOwnedSlice(),
                    .fys = fys.toOwnedSlice(),
                    .blocks = blocks.toOwnedSlice(),
                };
            }

            pub fn deinit(self: *Self, alloc: Allocator) void {
                alloc.free(self.fxs);
                alloc.free(self.fys);
                alloc.free(self.blocks);
            }
        };

        return Blocks;
    }
};

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

    pub fn to_string(self: BigCoord, alloc: Allocator) ![]u8 {
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

    /// iter counts (0 if not calculated, +ve = iter count, -ve = last failed count)
    iters: []i32,

    pub fn init(self: *Self, arena: Allocator, sz: u16) !void {
        self.sz = sz;
        self.iters = try arena.alloc(i32, @intCast(u32, sz) * sz);
        std.mem.set(i32, self.iters, 0);
    }

    pub fn deinit(self: *Self, arena: Allocator) void {
        arena.free(self.iters);
    }

    pub inline fn iter(self: Self, x: u32, y: u32) i32 {
        return self.iters[y * self.sz + x];
    }

    pub inline fn setIter(self: Self, x: u32, y: u32, it: i32) void {
        self.iters[y * self.sz + x] = it;
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
    pub fn get(self: Self, coord: BigCoord) ?*BlockData {
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
