const std = @import("std");

var random = std.rand.SplitMix64.init(1);

pub const XY = struct { x: u32, y: u32 };

fn walkSegments2(xys: []XY, SEGS: u32) void {
    const HSEGS = @divExact(SEGS, 2);
    var idx: u32 = 0;
    {
        // top
        var byi: u32 = 0;
        while (byi < HSEGS) : (byi += 1) {
            const by = HSEGS - byi - 1;
            {
                // left
                var bx: u32 = 0;
                while (bx < HSEGS) : (bx += 1) {
                    xys[idx] = .{ .x = HSEGS - bx - 1, .y = by };
                    idx += 1;
                }
            }
            {
                // right
                var bx: u32 = 0;
                while (bx < HSEGS) : (bx += 1) {
                    xys[idx] = .{ .x = HSEGS + bx, .y = by };
                    idx += 1;
                }
            }
        }
    }
    {
        // bottom
        var byi: u32 = 0;
        while (byi < HSEGS) : (byi += 1) {
            const by = HSEGS + byi;
            {
                // left
                var bx: u32 = 0;
                while (bx < HSEGS) : (bx += 1) {
                    xys[idx] = .{ .x = HSEGS - bx - 1, .y = by };
                    idx += 1;
                }
            }
            {
                // right
                var bx: u32 = 0;
                while (bx < HSEGS) : (bx += 1) {
                    xys[idx] = .{ .x = HSEGS + bx, .y = by };
                    idx += 1;
                }
            }
        }
    }
}

fn fillBoxesInnerToOuterRange(xys: []XY, idx0: u32, cx: u32, cy: u32, s: u32) u32 {
    if (s == 0) {
        return idx0;
    }

    // draw smaller boxes
    var idx = idx0;
    if (s > 0) {
        idx = fillBoxesInnerToOuterRange(xys, idx, cx, cy, s - 2);
    }

    // std.debug.print("s={}\n", .{s});
    const hs = @divExact(s, 2);
    var iy: u32 = 0;
    while (iy < s) : (iy += 1) {
        var ix: u32 = 0;
        while (ix < s) : (ix += 1) {
            if ((ix != 0 and ix != s - 1) and (iy != 0 and iy != s - 1)) continue;
            // const x = cx + ix - hs;
            const x = if (iy == 0 or (iy & 1) == 0) cx + ix - hs else cx + (s - 1 - ix) - hs;
            const y = cy + iy - hs;
            // std.debug.print("  {},{}\n", .{x,y});
            xys[idx] = .{ .x = x, .y = y };
            idx += 1;
        }
    }

    return idx;
}

pub fn fillBoxesInnerToOuter(xys: []XY, SEGS: u32) void {
    var idx: u32 = 0;
    idx = fillBoxesInnerToOuterRange(xys, idx, SEGS / 2, SEGS / 2, SEGS);
    if (idx != SEGS * SEGS) unreachable;
}

fn walkSegmentsTL(xys: []XY, SEGS: u32) void {
    var idx: u32 = 0;
    var by: u32 = 0;
    while (by < SEGS) : (by += 1) {
        var bx: u32 = 0;
        while (bx < SEGS) : (bx += 1) {
            xys[idx] = .{ .x = bx, .y = by };
            idx += 1;
        }
    }
}

/// Shuffle the boxes for less predictable rendering
pub fn shuffle(xys: []XY) void {
    const NSEGS = xys.len;
    var idx : u32 = 0;
    var xy: XY = undefined;
    while (idx < NSEGS / 2) {
        const dest = random.next() % NSEGS;
        if (idx == dest) continue;
        // std.debug.print("{} <-> {}\n", .{ idx, dest });
        xy = xys[idx];
        xys[idx] = xys[dest];
        xys[dest] = xy;
        idx += 1;
    }

}
