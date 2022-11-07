const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;
const json = std.json;

const bignum = @import("bignum.zig");
const Mandel = bignum.Mandel;

const boxfill = @import("boxfill.zig");
const XY = boxfill.XY;

const rb = @import("rb.zig");

pub const Params = struct {
    sx: u32,
    sy: u32,
    blockSize: u16,
    zoom: u16,
    words: u16,
    magShift: u5,
    iters: u32,
    cx: []u8,
    cy: []u8,

    pub fn span(self: Params) bignum.NativeFloat {
        return std.math.scalbn(@floatCast(bignum.NativeFloat, 4.0), -@intCast(i16, self.zoom));
    }

    pub fn getDefaultIntSize(self: Params) u16 {
        const zoom = self.zoom;
        const size: u16 = (zoom + 63) >> 6;
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

    pub fn readConfig(self: *Params, alloc: Allocator, path: []const u8) !void {
        var file = try std.fs.cwd().openFile(path, .{});
        var parser = std.json.Parser.init(alloc, false);
        defer parser.deinit();
        // try json.stringify(self, .{}, file.reader());
        var input = try alloc.alloc(u8, 1024);
        defer alloc.free(input);
        var size = try file.reader().readAll(input);
        var text = input[0..size];
        std.debug.print("text = {s}\n", .{text});

        var tree = try parser.parse(text);
        // self.*.cx = try parseHexFloat(bignum.NativeFloat, tree.root.Object.get("cx").?.String);
        // self.*.cy = try parseHexFloat(bignum.NativeFloat, tree.root.Object.get("cy").?.String);
        alloc.free(self.*.cx);
        alloc.free(self.*.cy);
        self.*.cx = try alloc.dupe(u8, tree.root.Object.get("cx").?.String);
        self.*.cy = try alloc.dupe(u8, tree.root.Object.get("cy").?.String);
        self.*.zoom = @intCast(u16, tree.root.Object.get("zoom").?.Integer);
        self.*.iters = @intCast(u32, tree.root.Object.get("iters").?.Integer);
        self.*.magShift = @intCast(u5, tree.root.Object.get("magShift").?.Integer);
        const jsonBits = tree.root.Object.get("words");
        self.*.words = if (jsonBits) |j| @intCast(u16, j.Integer) else self.*.getDefaultIntSize();
        // self.*.sx = tree.root.Object.get("sx").?.Float;
        // self.*.sy = tree.root.Object.get("sy").?.Float;
    }

    pub fn BlockWorkMaker(comptime T: type) type {
        const BlockCalc = struct {
            const Self = @This();
            const BigFloat = bignum.BigFixedFloat(T, 8);

            data: *BlockData,

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

            pub fn calculate(self: Self, stop: *std.atomic.Atomic(bool)) !void {
                const rectSize = self.rectSize;
                const data = self.data;

                var all : u32 = 0;
                var recalc : u32 = 0;
                var oy: u32 = 0;
                while (!stop.*.load(.Unordered) and oy < data.sz) : (oy += rectSize) {
                    const py = self.py + oy;
                    const fy = self.fys[py];
                    var ox: u32 = 0;
                    while (ox < data.sz) : (ox += rectSize) {
                        const px = self.px + ox;

                        // already done?
                        const iterIndex = oy * data.sz + ox;
                        const current = data.iters[iterIndex];
                        if (current == 0 or (current < 0 and -current < self.maxIters)) {
                            // need to calculate (more)
                            const iters = BigFloat.calcMandelbrot(self.fxs[px], fy, self.maxIters);
                            data.iters[iterIndex] = iters;
                            recalc += 1;
                        }
                        all += 1;
                    }
                }

                if (all > 0) {
                    const perc = recalc * 16 / all;
                    std.debug.print("{s}", .{".123456789ABCDEF="[perc..perc+1]});
                }
            }
        };

        const Blocks = struct {
            const Self = @This();
            const BigFloatType = bignum.BigFixedFloat(T, 8);

            blocks: []BlockCalc,
            fxs: []T,
            fys: []T,

            pub fn init(alloc: Allocator, params: Params, layer: *MandelLayer) !Self {
                // var params: Params = selfOrig;

                var fxs = std.ArrayList(T).init(alloc);
                var fys = std.ArrayList(T).init(alloc);

                std.debug.print("span = {}\n", .{params.span()});

                const fx0 = try bignum.parseBig(T, params.cx);
                const fy0 = try bignum.parseBig(T, params.cy);

                // split drawscape into segments
                const ps = params.blockSize;
                const SEGS: u32 = @divExact(std.math.min(params.sx, params.sy), ps);

                std.debug.print("SEGS = {}, ps = {}\n", .{ SEGS, ps });

                // calculate X and Y coords once
                {
                    const fxStep = BigFloatType.fromFloat(std.math.scalbn(params.span(), -@intCast(i32, std.math.log2_int(u32, params.sx))));
                    const hx = @intCast(i32, params.sx / 2);
                    var px: i32 = -hx;
                    while (px < params.sx) : (px += 1) {
                        try fxs.append(fx0 + fxStep * px);
                    }
                }
                {
                    const fyStep = BigFloatType.fromFloat(std.math.scalbn(params.span(), -@intCast(i32, std.math.log2_int(u32, params.sy))));
                    const hy = @intCast(i32, params.sy / 2);
                    var py: i32 = -hy;
                    while (py < params.sy) : (py += 1) {
                        try fys.append(fy0 + fyStep * py);
                    }
                }


                // then make blocks that reference indices into those arrays
                var blocks = std.ArrayList(BlockCalc).init(alloc);

                const rectSize: u32 = @as(u32, 1) << params.magShift;

                // get a calculation order, for quickest results
                // in the interesting area of the center
                var xys: []XY = try alloc.alloc(XY, SEGS*SEGS);
                defer alloc.free(xys);
                boxfill.fillBoxesInnerToOuter(xys, SEGS);

                // per block...
                for (xys) |xy| {
                    const px0 = xy.x * ps;
                    const py0 = xy.y * ps;

                    const coord = try BigCoord.init(T, alloc, fxs.items[px0], fys.items[py0]);
                    const blockData = try layer.ensure(coord);

                    try blocks.append(.{
                        .px = px0,
                        .py = py0,
                        .maxIters = params.iters,
                        .rectSize = rectSize,
                        .fxs = fxs.items,
                        .fys = fys.items,
                        // .rects = try alloc.alloc(RectColor, pxs * pys),
                        .data = blockData,
                    });
                }

                std.debug.print("Made {} items\n", .{blocks.items.len});
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

/// Coordinate for marking the edge of a MandelViewpoint.
pub const BigCoord = struct {
    const Self = @This();

    x: []u64,
    y: []u64,

    pub fn init(comptime IntType: type, alloc: Allocator, xi : IntType, yi: IntType) !Self {
        return .{
            .x = try alloc.dupe(u64, bignum.intWords(IntType, xi)),
            .y = try alloc.dupe(u64, bignum.intWords(IntType, yi)),
        };
    }

    pub fn hash(self: Self) u64 {
        var hasher = std.hash.Wyhash.init(0);
        // .Deep means, go into the slice, rather than taking the address
        std.hash.autoHashStrat(&hasher, self.x, .Deep);
        std.hash.autoHashStrat(&hasher, self.y, .Deep);
        const h = hasher.final();
        return h;
    }
    pub fn eql(self: BigCoord, other: BigCoord) bool {
        return std.mem.eql(u64, self.x, other.x) and std.mem.eql(u64, self.y, other.y);
    }

    pub fn to_string(self: BigCoord, alloc: Allocator) ![]u8 {
        var list = std.ArrayList(u8).init(alloc);
        try list.appendSlice("[");
        var xs = try bignum.printBigArray(alloc, self.x); defer alloc.free(xs);
        try list.appendSlice(xs);
        try list.appendSlice(", ");
        var ys = try bignum.printBigArray(alloc, self.y); defer alloc.free(ys);
        try list.appendSlice(ys);
        try list.appendSlice("]");
        return list.toOwnedSlice();
    }
};

/// Calculated block data, meaningful relative to some Params
pub const BlockData = struct {
    const Self = @This();

    /// coord of box UL
    coord: BigCoord,

    /// size of box in pixels
    sz: u16,

    /// iter counts (0 if not calculated, +ve = iter count, -ve = last failed count)
    iters: []i32,

    pub fn init(self: *Self, alloc: Allocator, coord: BigCoord, sz: u16) !void {
        self.coord = coord;
        self.sz = sz;
        self.iters = try alloc.alloc(i32, sz * sz);
        std.mem.set(i32, self.iters, 0);
    }

    pub fn deinit(self: @This(), alloc: Allocator) void {
        alloc.free(self.iters);
        // coord not owned
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

    alloc: Allocator,
    zoom: u16,
    /// block size
    sz: u16,
    blocks: BlockMap,

    pub fn init(self: *Self, alloc: Allocator, zoom: u16, sz: u16) !void {
        self.alloc = alloc;
        self.zoom = zoom;
        self.sz = sz;
        self.blocks = BlockMap.init(alloc);
    }

    pub fn deinit(self: *Self) void {
        self.blocks.clearAndFree();
        self.blocks.deinit();
    }

    /// Get the block
    pub fn ensure(self: *Self, coord: BigCoord) !*BlockData {
        var res = try self.blocks.getOrPut(coord);
        var blockPtr = res.value_ptr;
        if (!res.found_existing) {
            blockPtr.* = try self.alloc.create(BlockData);
            try blockPtr.*.init(self.alloc, coord, self.sz);
        }
        return blockPtr.*;
    }
};

/// All the data we've computed
pub const MandelStorage = struct {
    const Self = @This();
    const LayerMap = std.HashMap(u16, MandelLayer, struct {
        pub fn hash(_: anytype, a: u16) u64 { return a; }
        pub fn eql(_: anytype, a: u16, b: u16) bool { return a == b; }
    }, 80);

    alloc: Allocator,
    layers: LayerMap,

    pub fn init(alloc: Allocator) !Self {
        return .{
            .alloc = alloc,
            .layers = LayerMap.init(alloc)
        };
    }

    pub fn deinit(self: *Self) void {
        self.layers.clearAndFree();
        self.layers.deinit();
    }

    pub fn ensure(self: *Self, zoom: u16, sz: u16) !*MandelLayer {
        var gop = try self.layers.getOrPut(zoom);
        var layer = gop.value_ptr;
        if (!gop.found_existing) {
            std.debug.print("new zoom level layer {}\n", .{zoom});
            try layer.init(self.alloc, zoom, sz);
        }
        return layer;
    }
};
