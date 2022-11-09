const std = @import("std");
const viewer = @import("viewer.zig");

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{
        // .verbose_log=true
    }){};
    var alloc = gpa.allocator();
    // const alloc = std.heap.c_allocator;
    var v = try viewer.Viewer.init(alloc, 1024, 2048, 128);
    try v.run();
    std.debug.print("Bye!\n", .{});
    v.deinit();
}
