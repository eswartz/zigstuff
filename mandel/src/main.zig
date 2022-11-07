const std = @import("std");
const viewer = @import("viewer.zig");

pub fn main() !void {
    var v = try viewer.Viewer.init(std.heap.page_allocator);
    try v.run();
    std.debug.print("Bye!\n", .{});
    v.deinit();
}
