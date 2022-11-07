const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;
const json = std.json;

const bignum = @import("bignum.zig");
const Mandel = bignum.Mandel;

const boxfill = @import("boxfill.zig");
const XY = boxfill.XY;

const sdl2 = @cImport({
    @cInclude("SDL2/SDL.h");
});

const RectColor = struct {
    rect: sdl2.SDL_FRect,
    cnt: u32,
};

const Params = struct {
    sx: u32,
    sy: u32,
    zoom: u16,
    words: u16,
    magShift: u5,
    iters: i32,
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

    fn writeConfig(self: Params, alloc: Allocator, path: []const u8) !void {
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

    fn readConfig(self: *Params, alloc: Allocator, path: []const u8) !void {
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
        self.*.iters = @intCast(i32, tree.root.Object.get("iters").?.Integer);
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

            params: Params,
            px: u32,
            py: u32,
            sx: u32,
            sy: u32,
            rectSize: u32,
            fxs: []T,
            fys: []T,
            rects: []RectColor,

            pub fn calculate(self: Self, stop: *std.atomic.Atomic(bool)) !void {
                const rectSize = self.rectSize;

                var idx: u32 = 0;
                var oy: u32 = 0;
                while (!stop.*.load(.Unordered) and oy < self.sy) : (oy += rectSize) {
                    const py = self.py + oy;
                    const fy = self.fys[py];
                    var ox: u32 = 0;
                    while (ox < self.sx) : (ox += rectSize) {
                        const px = self.px + ox;
                        const iters = BigFloat.calcMandelbrot(self.fxs[px], fy, self.params.iters);

                        var rectColor = &self.rects[idx];
                        idx += 1;
                        rectColor.*.rect.x = @intToFloat(f32, px);
                        rectColor.*.rect.y = @intToFloat(f32, py);
                        rectColor.*.rect.w = @intToFloat(f32, rectSize);
                        rectColor.*.rect.h = @intToFloat(f32, rectSize);
                        rectColor.*.cnt = @intCast(u32, if (iters >= 0) iters else 0);
                    }
                }
            }

            pub fn render(self: Self, r: ?*sdl2.SDL_Renderer) void {
                for (self.rects) |rect| {
                    var g: u8 = @intCast(u8, rect.cnt & 0xff);
                    g = (g << 7) | (g >> 1);
                    _ = sdl2.SDL_SetRenderDrawColor(r, g, g, g, 255);
                    _ = sdl2.SDL_RenderFillRectF(r, &rect.rect);
                }

                _ = sdl2.SDL_RenderPresent(r);
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
                    try blocks.append(.{
                        .params = self, // copy
                        .px = px0,
                        .py = py0,
                        .sx = pxs,
                        .sy = pys,
                        .rectSize = rectSize,
                        .fxs = fxs.items,
                        .fys = fys.items,
                        .rects = try alloc.alloc(RectColor, pxs * pys),
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
                    alloc.free(block.rects);
                }
                alloc.free(self.blocks);
            }
        };

        return Blocks;
    }
};

pub const MandelError = error{
    SDLError,
};

pub const Viewer = struct {
    const Self = @This();

    alloc: Allocator,
    window: *sdl2.SDL_Window,
    renderer: *sdl2.SDL_Renderer,
    params: Params,
    exit: bool,
    saveFile : []u8,

    pub fn init(alloc: Allocator) !Self {
        sdl2.SDL_SetMainReady();

        if (sdl2.SDL_Init(sdl2.SDL_INIT_VIDEO) != 0) return error.SDLError;
        errdefer sdl2.SDL_Quit();
        if (sdl2.SDL_VideoInit(0) != 0) return error.SDLError;

        var window = sdl2.SDL_CreateWindow("Mandelbrot!", sdl2.SDL_WINDOWPOS_UNDEFINED, sdl2.SDL_WINDOWPOS_UNDEFINED, 1024, 1024, sdl2.SDL_WINDOW_RESIZABLE | sdl2.SDL_WINDOW_ALLOW_HIGHDPI);
        if (window == null) return error.SDLError;

        var renderer = sdl2.SDL_CreateRenderer(window, -1, 0);
        if (renderer == null) return error.SDLError;

        var params = Params{
            .sx = 1024,
            .sy = 1024,
            // .cx = 0.285, .cy = 0.01
            // .cx = 2.86953125e-01, .cy = 1.20751953125e-02,
            // .cx = 2.869521827343851e-01, .cy = 1.207519461400807e-02,
            // .cx = 2.8695218267617745e-01, .cy = 1.207519461400807e-02,
            // .cx = 2.86952182676179e-01, .cy = 1.2075194614007514e-02,
            // .cx = 2.86952182676179e-01, .cy = 1.20751946140075e-02,
            // .cx = try parseHexFloat(F, "0x1.25d6cb0070a58fb5cbb8a349469ep-2"),
            // .cy = try parseHexFloat(F, "0x1.8bae12fae1332fedf28b4701e9a0p-7"),
            // .cx = try parseHexFloat(F, "0x1.8p-0"),
            // .cy = try parseHexFloat(F, "0x1.aaae12fae1p-0"),
            // .zoom = 106, .magShift = 2,
            // .cx = "0x1.25d6cb0070370fb5cbb8a349469ep-2",
            // .cy = "0x1.8bae12fafb2b2fedf28b4701e9a0p-7",
            .cx = try alloc.dupe(u8, ".0000000000000000"),
            .cy = try alloc.dupe(u8, ".0000000000000000"),
            .zoom = 0,
            .magShift = 2,
            .words = 1,
            .iters = 300,
            // .iters = 400,
        };

        var saveName = try alloc.dupe(u8, "save.json");
        {
            var argIter = std.process.args();
            _ = argIter.next();
            if (argIter.next()) |arg| { saveName = try alloc.dupe(u8, arg); std.debug.print("save name = {s}\n", .{saveName}); }
        }

        return Viewer{
            .alloc = alloc,
            .window = window.?,
            .renderer = renderer.?,
            .params = params,
            .exit = false,
            .saveFile = saveName,
        };
    }

    pub fn run(self: *Self) !void {
        while (!self.exit) {
            try switch (self.params.words) {
                inline 1...bignum.MAXWORDS => |v| self.renderForParams(bignum.BigInt(64 * v)),
                else => self.renderForParams(bignum.BigInt(64 * bignum.MAXWORDS)),
            };
        }
    }

    pub fn deinit(self: Self) void {
        sdl2.SDL_DestroyRenderer(self.renderer);
        sdl2.SDL_DestroyWindow(self.window);
    }

    fn renderForParams(self: *Self, comptime BigIntType: type) !void {
        std.debug.print("Rendering with BigInt({s}) at x= {s}, y= {s}, zoom= {}, iters= {}...\n", .{ @typeName(BigIntType), self.params.cx, self.params.cy, self.params.zoom, self.params.iters });

        var timer = try std.time.Timer.start();
        const timeStart = timer.lap();

        var workLoad = try Params.BlockMaker(BigIntType).init(self.alloc, self.params);
        defer workLoad.deinit(self.alloc);

        std.debug.print("time for making blocks = {}\n", .{timer.lap() - timeStart});

        const BlockType = @TypeOf(workLoad.blocks[0]);

        const ThreadInfo = struct {
            fn calcThread(stop: *std.atomic.Atomic(bool), readyBlocks: *std.ArrayList(BlockType), readyMutex: *std.Thread.Mutex, blocksLeft: *std.atomic.Atomic(usize), blocks: []BlockType) !void {
                for (blocks) |block| {
                    try block.calculate(stop);

                    readyMutex.*.lock();
                    try readyBlocks.*.append(block);
                    if (!stop.*.load(.Unordered)) {
                        _ = blocksLeft.*.fetchSub(1, .Acquire);
                    }
                    readyMutex.*.unlock();
                }
            }
        };

        // make worker threads
        var stop = std.atomic.Atomic(bool).init(false);
        var readyBlocks = std.ArrayList(BlockType).init(self.alloc);
        defer readyBlocks.deinit();
        var readyMutex: std.Thread.Mutex = undefined;
        var blocksLeft = std.atomic.Atomic(usize).init(workLoad.blocks.len);

        const NTHREADS: usize = std.math.min(256, try std.Thread.getCpuCount());
        var threads: [256]std.Thread = undefined;

        const PERTHREAD = @divExact(workLoad.blocks.len, NTHREADS);

        var nt: u32 = 0;
        while (nt < NTHREADS) : (nt += 1) {
            const slice = workLoad.blocks[nt * PERTHREAD .. (nt + 1) * PERTHREAD];
            threads[nt] = try std.Thread.spawn(.{}, ThreadInfo.calcThread, .{ &stop, &readyBlocks, &readyMutex, &blocksLeft, slice });
        }

        std.debug.print("time for thread launch = {}\n", .{timer.lap() - timeStart});

        // render as work arrives
        var cont: bool = false;
        var done: bool = false;

        _ = sdl2.SDL_RenderClear(self.renderer);
        wait: while (!done) {
            cont = try self.handleInput();

            if (cont) {
                break;
            }

            readyMutex.lock();
            var block: ?BlockType = null;
            while (true) {
                block = readyBlocks.popOrNull();
                if (block == null) {
                    if (blocksLeft.load(.Acquire) == 0) {
                        done = true;
                        std.debug.print("time for calculation   = {d}s\n", .{@intToFloat(f64, (timer.lap() - timeStart)) / 1.0e9});
                        readyMutex.unlock();
                        break :wait;
                    }
                    break;
                }
                block.?.render(self.renderer);
            }
            readyMutex.unlock();
        }

        if (!self.exit) {
            while (!cont) {
                cont = try self.handleInput();
            }
        }

        stop.store(true, .Unordered);
        for (threads[0..NTHREADS]) |thread| {
            thread.join();
        }
    }

    fn handleInput(self: *Self) !bool {
        var ps = &self.params;
        var event: sdl2.SDL_Event = undefined;
        var cont = false;
        while (sdl2.SDL_WaitEventTimeout(&event, 10) != 0) {
            var words = ps.*.words;
            var span = ps.*.span();
            if (event.type == sdl2.SDL_QUIT) {
                self.exit = true;
            } else if (event.type == sdl2.SDL_KEYDOWN) {
                var ke = @ptrCast(*sdl2.SDL_KeyboardEvent, &event);
                cont = true;
                const shiftAmt: u6 = if ((ke.keysym.mod & sdl2.KMOD_LSHIFT) != 0) 4 else 1;
                _ = switch (ke.keysym.sym) {
                    sdl2.SDLK_ESCAPE => self.exit = true,
                    sdl2.SDLK_PLUS, sdl2.SDLK_EQUALS => ps.*.zoom += 1,
                    sdl2.SDLK_MINUS => ps.*.zoom = if (ps.*.zoom > 1) ps.*.zoom - 1 else 0,
                    sdl2.SDLK_UP => try bignum.faddShift(self.alloc, &ps.*.cy, words, -span, shiftAmt),
                    sdl2.SDLK_DOWN => try bignum.faddShift(self.alloc, &ps.*.cy, words, span, shiftAmt),
                    sdl2.SDLK_LEFT => try bignum.faddShift(self.alloc, &ps.*.cx, words, -span, shiftAmt),
                    sdl2.SDLK_RIGHT => try bignum.faddShift(self.alloc, &ps.*.cx, words, span, shiftAmt),
                    sdl2.SDLK_PAGEDOWN => ps.*.iters = std.math.max(1, ps.*.iters -% 50),
                    sdl2.SDLK_PAGEUP => ps.*.iters += 50,
                    sdl2.SDLK_LEFTBRACKET => ps.*.words = if (words > 1) words - 1 else 1,
                    sdl2.SDLK_RIGHTBRACKET => ps.*.words += 1,
                    sdl2.SDLK_SPACE => {
                        try bignum.falign(self.alloc, &ps.*.cx, words);
                        try bignum.falign(self.alloc, &ps.*.cy, words);
                    },
                    sdl2.SDLK_F5 => {
                        try ps.*.writeConfig(self.alloc, self.saveFile);
                        cont = false;
                    },
                    sdl2.SDLK_F9 => {
                        if (ps.*.readConfig(self.alloc, self.saveFile)) {} else |err| {
                            std.debug.print("Failed to load: {}\n", .{err});
                            cont = false;
                        }
                    },
                    '0'...'9' => |v| ps.*.magShift = @intCast(u5, v - '0'),
                    else => cont = false,
                };
            } else if (event.type == sdl2.SDL_MOUSEBUTTONDOWN) {
                var me = @ptrCast(*sdl2.SDL_MouseButtonEvent, &event);
                if (me.button == 1) {
                    cont = true;
                    // recenter on mouse position
                    try bignum.fadj(self.alloc, &ps.*.cy, words, span, me.y - @intCast(i32, ps.*.sy >> 1), ps.*.sy);
                    try bignum.fadj(self.alloc, &ps.*.cx, words, span, me.x - @intCast(i32, ps.*.sx >> 1), ps.*.sx);
                }
            }
        }
        return cont;
    }
};
