const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;
const json = std.json;

const bignum = @import("bignum.zig");
const Mandel = bignum.Mandel;

const boxfill = @import("boxfill.zig");
const XY = boxfill.XY;

pub const Params = struct {
    sx: u32,
    sy: u32,
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

    pub fn BlockMaker(comptime T: type) type {
        const Block = struct {
            const Self = @This();
            const BigFloat = bignum.BigFixedFloat(T, 8);

            maxIters: u32,
            px: u32,
            py: u32,
            sx: u32,
            sy: u32,
            rectSize: u32,
            // reference copy
            fxs: []T,
            // reference copy
            fys: []T,

            data: BlockData,

            pub fn calculate(self: Self, stop: *std.atomic.Atomic(bool)) !void {
                const rectSize = self.rectSize;

                // var idx: u32 = 0;
                var oy: u32 = 0;
                while (!stop.*.load(.Unordered) and oy < self.sy) : (oy += rectSize) {
                    const py = self.py + oy;
                    const fy = self.fys[py];
                    var ox: u32 = 0;
                    while (ox < self.sx) : (ox += rectSize) {
                        const px = self.px + ox;

                        // already done?
                        const iterIndex = oy * self.sx + ox;
                        const current = self.data.iters[iterIndex];
                        if (current == 0 or (current < 0 and -current < self.maxIters)) {
                            // need to calculate (more)
                            const iters = BigFloat.calcMandelbrot(self.fxs[px], fy, self.maxIters);
                            self.data.iters[iterIndex] = iters;
                        }

                        // var rectColor = &self.rects[idx];
                        // idx += 1;
                        // rectColor.*.rect.x = @intToFloat(f32, px);
                        // rectColor.*.rect.y = @intToFloat(f32, py);
                        // rectColor.*.rect.w = @intToFloat(f32, rectSize);
                        // rectColor.*.rect.h = @intToFloat(f32, rectSize);
                        // rectColor.*.cnt = @intCast(u32, if (iters >= 0) iters else 0);
                    }
                }
            }

            pub fn deinit(self: Self, alloc: Allocator) void {
                self.data.deinit(alloc);
            }
        };

        const Blocks = struct {
            const Self = @This();
            const BigFloatType = bignum.BigFixedFloat(T, 8);

            blocks: []Block,
            fxs: []T,
            fys: []T,

            pub fn init(alloc: Allocator, selfOrig: Params) !Self {
                var self: Params = selfOrig;

                var fxs = std.ArrayList(T).init(alloc);
                var fys = std.ArrayList(T).init(alloc);

                std.debug.print("span = {}\n", .{self.span()});
                // calculate X and Y coords once
                {
                    const fx0 = try bignum.parseBig(T, self.cx);
                    const fxStep = BigFloatType.fromFloat(std.math.scalbn(self.span(), -@intCast(i32, std.math.log2_int(u32, self.sx))));
                    const hx = @intCast(i32, self.sx / 2);
                    var px: i32 = -hx;
                    while (px < self.sx) : (px += 1) {
                        try fxs.append(fx0 + fxStep * px);
                    }
                }
                {
                    const fy0 = try bignum.parseBig(T, self.cy);
                    const fyStep = BigFloatType.fromFloat(std.math.scalbn(self.span(), -@intCast(i32, std.math.log2_int(u32, self.sy))));
                    const hy = @intCast(i32, self.sy / 2);
                    var py: i32 = -hy;
                    while (py < self.sy) : (py += 1) {
                        try fys.append(fy0 + fyStep * py);
                    }
                }

                // then make blocks that reference indices into those arrays
                var blocks = std.ArrayList(Block).init(alloc);

                const rectSize: u32 = @as(u32, 1) << self.magShift;

                // split drawscape into segments
                const SEGS: u32 = 8;
                const pxs = @divExact(self.sx, SEGS);
                const pys = @divExact(self.sy, SEGS);

                var xys: [SEGS * SEGS]XY = undefined;
                boxfill.fillBoxesInnerToOuter(&xys, SEGS);

                for (xys) |xy| {
                    const px0 = xy.x * pxs;
                    const py0 = xy.y * pys;
                    var blockData = BlockData{
                        .px = px0,
                        .py = py0,
                        .sx = pxs,
                        .sy = pys,
                        .iters = undefined
                    };
                    try blockData.init(alloc);
                    try blocks.append(.{
                        .maxIters = self.iters,
                        .px = px0,
                        .py = py0,
                        .sx = pxs,
                        .sy = pys,
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
                for (self.blocks) |block| {
                    block.deinit(alloc);
                }
                alloc.free(self.blocks);
            }
        };

        return Blocks;
    }
};

/// Calculated block data, meaningful relative to some Params
pub const BlockData = struct {
    /// pixel coord of box (relative to params.{cx,cy})
    px: u32,
    py: u32,

    /// size of box (absolute, ignoring magnification)
    sx: u32,
    sy: u32,

    /// iter counts (for entire span of possible pixels;
    /// some may be unset (0) if rendered under magnification)
    iters: []i32,

    pub fn init(self: *@This(), alloc: Allocator) !void {
        self.iters = try alloc.alloc(i32, self.sx * self.sy);
        std.mem.set(i32, self.iters, 0);
    }
    pub fn deinit(self: @This(), alloc: Allocator) void {
        alloc.free(self.iters);
    }
};
