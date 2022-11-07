const std = @import("std");
const mandel = @import("mandel.zig");

pub fn main() !void {
    var viewer = try mandel.Viewer.init(std.heap.page_allocator);
    try viewer.run();
    viewer.deinit();
}
