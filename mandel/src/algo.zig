const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;
const json = std.json;

const config = @import("config.zig");
const BigIntMax = config.BigIntMax;

const bignum = @import("bignum.zig");
const BigInt = bignum.BigInt;
const BigFixedFloat = bignum.BigFixedFloat;
const Mandel = bignum.Mandel;

const cache = @import("cache.zig");
const BigCoord = cache.BigCoord;
const BlockData = cache.BlockData;
const MandelStorage = cache.MandelStorage;

const boxfill = @import("boxfill.zig");

const XY = @import("types.zig").XY;

pub const MandelAlgo = struct {
    pub fn init(comptime BigIntType: type) type {
        return struct {
            const BigFloat = BigFixedFloat(BigIntType, 8);
            const BigFloatType = BigFloat.Repr;
            const BigFloat2 = BigInt(@typeInfo(BigIntType).Int.bits * 2);

            pub fn mul(x: BigIntType, y: BigIntType) BigIntType {
                const ps = (@as(BigFloat2, x) * @as(BigFloat2, y)) >> BigFloat.mantissaBits;
                return @intCast(BigIntType, ps);
            }

            /// Calculate the intersection of a coordinate (creal,cimag) in the
            /// Mandelbrot set.  If it is not in the set, return a positive
            /// integer indicating the iteration count before exit.
            /// Else, return a negative number indicating the negative of how
            /// many iterations were attempted.
            pub fn calcMandelbrot(creal: BigFloatType, cimag: BigFloatType, maxIters: u32) i32 {
                const four: BigFloatType = BigFloat.fromFloat(4.0);
                var iters: i32 = @intCast(i32, maxIters);
                var real = creal;
                var imag = cimag;

                while (iters > 0) : (iters -= 1) {
                    const realsq = mul(real, real);
                    const imagsq = mul(imag, imag);
                    if (realsq + imagsq >= four) {
                        return @intCast(i32, maxIters) - iters + 1;
                    }

                    const temp = mul(real, imag);

                    // r = (r'*r' - i*i) + cr
                    real = (realsq - imagsq) +% creal;
                    // i = (r'*i' * 2) + ci
                    imag = (temp + temp) +% cimag;
                }
                return -@intCast(i32, maxIters);
            }
        };
    }
};

test "mb256_0" {
    const BigInt256 = BigInt(256);
    const m = BigFixedFloat(BigInt256, 8);
    const algo = MandelAlgo.init(BigInt256);
    {
        const x0 = m.fromFloat(-1.5);
        const y0 = m.fromFloat(1.5);
        try testing.expectEqual(@as(i32, 1), algo.calcMandelbrot(x0, y0, 100));
    }
    {
        const x0 = m.fromFloat(-1.25);
        const y0 = m.fromFloat(0.2);
        try testing.expectEqual(@as(i32, 11), algo.calcMandelbrot(x0, y0, 100));
    }

    const cx = -1.25;
    const cy = 0.2;

    const XS = 80;
    const YS = 32;
    const xrange = 0.1 / (XS + 0.0);
    const yrange = 0.1 / (YS + 0.0);

    var itersTot: i64 = 0;
    std.debug.print("\n", .{});
    var py: u32 = 0;
    while (py < YS) : (py += 1) {
        var px: u32 = 0;
        while (px < XS) : (px += 1) {
            const fx = m.fromFloat(cx + @intToFloat(bignum.NativeFloat, px) * xrange);
            const fy = m.fromFloat(cy + @intToFloat(bignum.NativeFloat, py) * yrange);
            const iters = algo.calcMandelbrot(fx, fy, 256);
            itersTot += if (iters > 0) iters else -iters;
            std.debug.print("{c}", .{if (iters >= 0) "-ABCDEFGHIJKLMNOPQRSTUVWXYZ"[@intCast(u32, iters) % 26] else ' '});
        }
        std.debug.print("\n", .{});
    }
    try testing.expectEqual(@as(i64, 98465), itersTot);
}

pub const Params = struct {
    const ZoomBitMap = std.AutoArrayHashMap(u16, u16);
    const ZoomIterMap = std.AutoArrayHashMap(u16, u32);

    sx: u32,
    sy: u32,
    zoom: u16,
    words: u16,
    magShift: u5,
    iters: u32,
    cx: []u8,
    cy: []u8,
    zoomBits: ZoomBitMap,
    zoomIters: ZoomIterMap,

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
            .zoomBits = ZoomBitMap.init(alloc),
            .zoomIters = ZoomIterMap.init(alloc),
        };
    }

    pub fn deinit(self: *Params, alloc: Allocator) void {
        alloc.free(self.cx);
        alloc.free(self.cy);
        self.zoomBits.deinit();
        self.zoomIters.deinit();
    }

    pub fn setCx(self: *Params, alloc: Allocator, str: []u8) void {
        alloc.free(self.cx);
        self.cx = str;
    }
    pub fn setCy(self: *Params, alloc: Allocator, str: []u8) void {
        alloc.free(self.cy);
        self.cy = str;
    }

    pub fn span(self: Params, calcSize: u32) bignum.NativeFloat {
        const fullSpan = std.math.scalbn(@floatCast(bignum.NativeFloat, 4.0), -@intCast(i16, self.zoom));
        return if (calcSize <= self.sx)
            fullSpan / @intToFloat(bignum.NativeFloat, @divExact(self.sx, calcSize))
        else
            fullSpan * @intToFloat(bignum.NativeFloat, @divExact(calcSize, self.sx));
    }

    pub fn segments(_: Params, calcSize: u32, blockSize: u16) u32 {
        return @divExact(calcSize, blockSize);
    }

    pub fn setZoom(self: *Params, zoom: u16) void {
        self.zoom = zoom;
        self.words = self.getIntSizeForZoom();
        self.iters = self.getItersForZoom();
    }

    pub fn setZoomInteractive(self: *Params, zoom: u16) void {
        const curZoom = self.zoom;
        const curWords = self.words;
        const curIters = self.iters;
        self.setZoom(zoom);
        if (curZoom < zoom) {
            // see if we fell off the end of the map and restore larger values
            if (self.words < curWords) {
                // probably fell off the end
                self.words = curWords;
            }
            if (self.iters < curIters) {
                self.iters = curIters;
            }
        }
    }

    /// Get number of words needed
    pub fn getDefaultIntSize(self: Params) u16 {
        const zoom = self.zoom;
        const size: u16 = 1 + ((zoom + std.math.log2_int(u32, self.sx)) >> 6);
        std.debug.print("zoom={}, isize={}\n", .{ zoom, size });
        return size;
    }

    /// Get number of words configured for zoom
    pub fn getIntSizeForZoom(self: Params) u16 {
        if (self.zoomBits.count() != 0) {
            const zoom = self.zoom;
            var it = self.zoomBits.iterator();
            var bits: u16 = 64;
            while (it.next()) |ent| {
                if (ent.key_ptr.* > zoom) break;
                bits = ent.value_ptr.*;
            }
            return (bits + 63) >> 6;
        }
        return self.getDefaultIntSize();
    }

    /// Remember current zoom and bits in zoomBits
    pub fn recordZoomBits(self: *Params) !void {
        std.debug.print("Remembering zoom bits {} for level {}\n", .{self.words << 6, self.zoom});
        try self.zoomBits.put(self.zoom, self.words << 6);
        self.sortZoomBits();
    }

    pub fn sortZoomBits(self: *Params) void {
        const Sorter = struct {
            entries: *ZoomBitMap.DataList,
            pub fn lessThan(s: @This(), ai: usize, bi: usize) bool {
                return s.entries.get(ai).key < s.entries.get(bi).key;
            }
        };
        self.zoomBits.sort(Sorter{ .entries = &self.zoomBits.unmanaged.entries });
    }

    /// Get iters configured for zoom
    pub fn getItersForZoom(self: Params) u32 {
        if (self.zoomIters.count() != 0) {
            const zoom = self.zoom;
            var it = self.zoomIters.iterator();
            var prevZoom: u16 = 0;
            var prevIters: u32 = 50 * zoom;
            var nextZoom: u16 = 0;
            var nextIters: u32 = 50 * zoom;
            while (it.next()) |ent| {
                prevZoom = nextZoom;
                prevIters = nextIters;
                nextZoom = ent.key_ptr.*;
                nextIters = ent.value_ptr.*;
                if (nextZoom > zoom) break;
            }
            return prevIters + (self.zoom - prevZoom) * (nextIters - prevIters) / (nextZoom - prevZoom);
        }
        return 50 * self.zoom;
    }

    /// Remember current zoom and iters in zoomIters
    pub fn recordZoomIters(self: *Params) !void {
        std.debug.print("Remembering zoom iters {} for level {}\n", .{self.iters, self.zoom});
        try self.zoomIters.put(self.zoom, self.iters);
        self.sortZoomIters();
    }

    pub fn sortZoomIters(self: *Params) void {
        const Sorter = struct {
            entries: *ZoomIterMap.DataList,
            pub fn lessThan(s: @This(), ai: usize, bi: usize) bool {
                return s.entries.get(ai).key < s.entries.get(bi).key;
            }
        };
        self.zoomIters.sort(Sorter{ .entries = &self.zoomIters.unmanaged.entries });
    }

    pub const BlockList = anyopaque;

    fn Block(comptime T: type) type {
        return struct {
            const Self = @This();
            const Algo = MandelAlgo.init(T);

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

            fn preFill(self: Self) void {
                const rectSize = self.rectSize;
                const data = self.data;

                const blockSize = data.sz;
                var oy: u32 = 0;
                const halfBlockSize = blockSize >> 1;

                // if we can, fill up from other levels first
                while (oy < blockSize) : (oy += rectSize) {
                    var ox: u32 = 0;
                    while (ox < blockSize) : (ox += rectSize) {
                        var iters = data.iter(ox, oy);
                        if (iters == 0) {
                            if (self.datam1 != null and (ox & 1) == 0 and (oy & 1) == 0) {
                                // reuse data from outer zoom
                                iters = self.datam1.?.iter(ox / 2, oy / 2);
                            }
                            if (iters == 0) {
                                // reuse data from inner zoom?
                                var inner: ?*const BlockData =
                                    if (ox < halfBlockSize)
                                    if (oy < halfBlockSize) self.datap1_00 else self.datap1_01
                                else if (oy < halfBlockSize) self.datap1_10 else self.datap1_11;
                                if (inner != null) {
                                    iters = inner.?.iter(((ox << 1) & (blockSize - 1)), ((oy << 1) & (blockSize - 1)));
                                }
                            }
                            data.setIter(ox, oy, iters);
                        }
                    }
                }
            }

            pub fn calculate(self: Self, stop: *std.atomic.Atomic(bool)) !bool {
                // if we can, fill up from other levels first
                if (self.datam1 != null or self.datap1_00 != null) {
                    self.preFill();
                }

                const rectSize = self.rectSize;
                const data = self.data;

                const blockSize = data.sz;
                var all: u32 = 0;
                var recalc: u32 = 0;
                var oy: u32 = 0;

                // inner loop
                while (!stop.*.load(.Unordered) and oy < blockSize) : (oy += rectSize) {
                    const py = self.py + oy;
                    const fy = self.fys[py];
                    var ox: u32 = 0;
                    while (ox < blockSize) : (ox += rectSize) {
                        const px = self.px + ox;

                        var iters = data.iter(ox, oy);
                        if (iters == 0 or (iters < 0 and -iters < self.maxIters)) {
                            // need to calculate (more)
                            iters = Algo.calcMandelbrot(self.fxs[px], fy, self.maxIters);
                            data.setIter(ox, oy, iters);
                            recalc += 1;
                        }
                        all += 1;
                    }
                }

                if (all > 0) {
                    const perc = recalc * 16 / all;
                    std.debug.print("{s}", .{".123456789ABCDEF="[perc .. perc + 1]});
                }
                return recalc > 0;
            }

            pub fn isComplete(self : Self) bool {
                const rectSize = self.rectSize;
                const data = self.data;

                const blockSize = data.sz;
                var oy: u32 = 0;

                while (oy < blockSize) : (oy += rectSize) {
                    var ox: u32 = 0;
                    while (ox < blockSize) : (ox += rectSize) {
                        const iters = data.iter(ox, oy);
                        if (iters == 0 or (iters < 0 and -iters < self.maxIters)) {
                            // need to calculate (more)
                            return false;
                        }
                    }
                }

                return true;
            }
        };
    }

    /// Create the blocks visible with the current MandelParams and
    /// a given maximal power-of-two window size.  This will add data to
    /// any currently-existing layers at the current zoom level and
    /// make use of the zoom level N-1 and N+1 to derive previously
    /// calculated data.
    pub fn BlockGenerator(comptime T: type) type {
        return struct {
            const BlockT = Block(T);
            const Self = @This();
            const BigFloatType = bignum.BigFixedFloat(T, 8);

            blocks: []BlockT,
            fxs: []T,
            fys: []T,

            // calcSize may be smaller or larger than params.sx/sy
            pub fn init(alloc: Allocator, params: Params, storage: *MandelStorage, calcSize: u32) !Self {
                var fxs = std.ArrayList(T).init(alloc);
                var fys = std.ArrayList(T).init(alloc);

                const pspan = params.span(calcSize);
                std.debug.print("span = {}\n", .{pspan});

                // split drawscape into segments
                const blockSize = storage.blockSize;
                const halfBlockSize = blockSize >> 1;
                const SEGS = params.segments(calcSize, blockSize);
                const div = @intToFloat(bignum.NativeFloat, blockSize * SEGS);
                var step = pspan / div;
                // std.debug.print("segs = {}, span = {}, div = {}, step = {}\n", .{ SEGS, pspan, div, step });
                // per-pixel step
                const floatStep = BigFloatType.fromFloat(step);

                // get raw coords, aligned to available precision
                const fx0 = try bignum.parseBig(T, params.cx) & -floatStep;
                const fy0 = try bignum.parseBig(T, params.cy) & -floatStep;

                std.debug.print("Center ==> {s}, {s}\n", .{ try bignum.printBig(alloc, T, fx0), try bignum.printBig(alloc, T, fy0) });

                var layer = try storage.ensure(params.zoom, blockSize);
                var layerM1 = if (params.zoom > 0) storage.get(params.zoom - 1) else null;
                var layerP1 = storage.get(params.zoom + 1);

                // var fss = try bignum.printBig(alloc, T, floatStep);
                // defer alloc.free(fss);
                // std.debug.print("fss = {s}\n", .{fss});

                // calculate X and Y coords once
                {
                    const hx = @intCast(i32, calcSize / 2);
                    var px: i32 = -hx;
                    while (px < calcSize) : (px += 1) {
                        try fxs.append(fx0 + floatStep * px);
                    }
                }
                {
                    const hy = @intCast(i32, calcSize / 2);
                    var py: i32 = -hy;
                    while (py < calcSize) : (py += 1) {
                        try fys.append(fy0 + floatStep * py);
                    }
                }

                // then make blocks that reference indices into those arrays
                var blocks = std.ArrayList(BlockT).init(alloc);

                const rectSize: u32 = @as(u32, 1) << params.magShift;

                // get a calculation order, for quickest results
                // in the interesting area of the center
                var xys: []XY = try alloc.alloc(XY, SEGS * SEGS);
                defer alloc.free(xys);
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

                // std.debug.print("Made {} blocks\n", .{blocks.items.len});
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
    }
};
