const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;
const json = std.json;

const bignum = @import("bignum.zig");
const Mandel = bignum.Mandel;

const boxfill = @import("boxfill.zig");
const XY = boxfill.XY;

const mandel = @import("mandel.zig");
const Params = mandel.Params;

const sdl2 = @cImport({
    @cInclude("SDL2/SDL.h");
});

const RectColor = struct {
    rect: sdl2.SDL_FRect,
    cnt: u32,
};

pub const Viewer = struct {
    const Self = @This();

    alloc: Allocator,
    window: *sdl2.SDL_Window,
    renderer: *sdl2.SDL_Renderer,
    params: Params,
    exit: bool,
    saveFile : []u8,
    storage : mandel.MandelStorage,

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

        var storage = try mandel.MandelStorage.init(alloc);

        return Viewer{
            .alloc = alloc,
            .window = window.?,
            .renderer = renderer.?,
            .params = params,
            .exit = false,
            .saveFile = saveName,
            .storage = storage
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

    fn renderBlock(self: Self, block: *mandel.BlockData) void {
        var rect : sdl2.SDL_FRect = undefined;

        const rectSize: u32 = @as(u32, 1) << self.params.magShift;
        var oy: u32 = 0;
        while (oy < block.sz) : (oy += rectSize) {
            const py = block.py + oy;
            var ox: u32 = 0;
            while (ox < block.sz) : (ox += rectSize) {
                const px = block.px + ox;

                const iterIndex = oy * block.sz + ox;
                const current = block.iters[iterIndex];

                var g: u8 = @intCast(u8, if (current < 0) 0 else current & 0xff);
                g = (g << 7) | (g >> 1);
                _ = sdl2.SDL_SetRenderDrawColor(self.renderer, g, g, g, 255);

                rect.x = @intToFloat(f32, px);
                rect.y = @intToFloat(f32, py);
                rect.w = @intToFloat(f32, rectSize);
                rect.h = @intToFloat(f32, rectSize);
                _ = sdl2.SDL_RenderFillRectF(self.renderer, &rect);
            }
        }
    }

    fn renderForParams(self: *Self, comptime BigIntType: type) !void {
        std.debug.print("Rendering with BigInt({s}) at x= {s}, y= {s}, zoom= {}, iters= {}...\n", .{ @typeName(BigIntType), self.params.cx, self.params.cy, self.params.zoom, self.params.iters });

        var timer = try std.time.Timer.start();
        const timeStart = timer.lap();

        var workLoad = try Params.BlockWorkMaker(BigIntType).init(self.alloc, self.params, &self.storage);
        defer workLoad.deinit(self.alloc);

        std.debug.print("time for making blocks = {}\n", .{timer.lap() - timeStart});

        for (workLoad.blocks) |block| {
            self.renderBlock(block.data);
        }

        std.debug.print("time for first block draw = {}\n", .{timer.lap() - timeStart});

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
            const start = nt * PERTHREAD;
            const slice = workLoad.blocks[start..start+PERTHREAD];
            threads[nt] = try std.Thread.spawn(.{}, ThreadInfo.calcThread, .{ &stop, &readyBlocks, &readyMutex, &blocksLeft, slice });
        }

        std.debug.print("time for thread launch = {}\n", .{timer.lap() - timeStart});

        // render as work arrives
        var cont: bool = false;
        var done: bool = false;

        while (!done) {
            cont = try self.handleInput();

            if (cont) {
                break;
            }

            readyMutex.lock();
            var block: ?BlockType = null;
            var any = false;
            while (true) {
                block = readyBlocks.popOrNull();
                if (block == null) {
                    if (blocksLeft.load(.Acquire) == 0) {
                        done = true;
                        std.debug.print("time for calculation   = {d}s\n", .{@intToFloat(f64, (timer.lap() - timeStart)) / 1.0e9});
                    }
                    break;
                }
                self.renderBlock(block.?.data);
                any = true;
            }
            readyMutex.unlock();

            if (any) {
                _ = sdl2.SDL_RenderPresent(self.renderer);
            }
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
                    sdl2.SDLK_ESCAPE => { self.exit = true; std.debug.print("Cancelling...\n", .{}); },
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
