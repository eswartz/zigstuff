const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;
const json = std.json;

const config = @import("config.zig");

const bignum = @import("bignum.zig");
const Mandel = bignum.Mandel;

const boxfill = @import("boxfill.zig");
const XY = boxfill.XY;

const mandel = @import("mandel.zig");
const Params = mandel.Params;

const cimports = @import("cimports.zig");
const sdl2 = cimports.sdl2;

const files = @import("files.zig");

pub const Viewer = struct {
    const Self = @This();

    alloc: Allocator,
    window: *sdl2.SDL_Window,
    renderer: *sdl2.SDL_Renderer,
    params: Params,
    exit: bool,
    pause: bool,
    saveFile : []u8,
    dataFile : []u8,
    storage : mandel.MandelStorage,

    pub fn init(alloc: Allocator) !Self {
        var winSize : c_int = 1024;
        var frameSize : u32 = 1024;
        var blockSize : u16 = 128;

        var saveName : []u8 = "";
        var dataName : []u8 = "";
        {
            var argIter = std.process.args();
            _ = argIter.next(); // name
            while (argIter.next()) |arg| {
                if (std.mem.eql(u8, arg, "-h")) {
                    std.debug.print(
                        \\ help: [options] .json .dat
                        \\
                        \\ -s pot: render size
                        \\ -r pot: window size
                        \\ -b pot: block size
                        \\
                        , .{});
                    std.process.exit(0);
                } if (std.mem.eql(u8, arg, "-s")) {
                    frameSize = try std.fmt.parseInt(u32, argIter.next().?, 10);
                    if (!std.math.isPowerOfTwo(frameSize)) {
                        return error.SizeMustBePowerOfTwo;
                    }
                } else if (std.mem.eql(u8, arg, "-r")) {
                    winSize = try std.fmt.parseInt(c_int, argIter.next().?, 10);
                    if (!std.math.isPowerOfTwo(winSize)) {
                        return error.WindowMustBePowerOfTwo;
                    }
                } else if (std.mem.eql(u8, arg, "-b")) {
                    blockSize = try std.fmt.parseInt(u16, argIter.next().?, 10);
                    if (!std.math.isPowerOfTwo(blockSize)) {
                        return error.BlockSizeMustBePowerOfTwo;
                    }
                } else if (arg[0] != '-') {
                    if (saveName.len == 0)
                        saveName = try alloc.dupe(u8, arg)
                    else if (dataName.len == 0)
                        dataName = try alloc.dupe(u8, arg)
                    else
                        return error.UnexpectedArgument;
                }
            }
        }

        if (blockSize > winSize) blockSize = @intCast(u16, winSize);
        if (winSize > frameSize) frameSize = @intCast(u32, winSize);

        if (saveName.len == 0) saveName = try alloc.dupe(u8, "save.json");
        if (dataName.len == 0) dataName = try alloc.dupe(u8, "data/zoom10_000.dat");

        var storage = try mandel.MandelStorage.init(alloc, blockSize);

        var params = try Params.init(alloc);

        params.sx = frameSize;
        params.sy = frameSize;

        sdl2.SDL_SetMainReady();

        if (sdl2.SDL_Init(sdl2.SDL_INIT_VIDEO) != 0) return error.SDLError;
        errdefer sdl2.SDL_Quit();
        if (sdl2.SDL_VideoInit(0) != 0) return error.SDLError;

        var window = sdl2.SDL_CreateWindow("Mandelbrot!", sdl2.SDL_WINDOWPOS_UNDEFINED, sdl2.SDL_WINDOWPOS_UNDEFINED,
            winSize, winSize,
            sdl2.SDL_WINDOW_RESIZABLE | sdl2.SDL_WINDOW_ALLOW_HIGHDPI);
        if (window == null) return error.SDLError;

        var renderer = sdl2.SDL_CreateRenderer(window, -1, 0);
        if (renderer == null) return error.SDLError;

        _ = sdl2.SDL_SetRenderDrawBlendMode(renderer, sdl2.SDL_BLENDMODE_BLEND);

        return Viewer{
            .alloc = alloc,
            .window = window.?,
            .renderer = renderer.?,
            .params = params,
            .exit = false,
            .pause = false,
            .saveFile = saveName,
            .dataFile = dataName,
            .storage = storage
        };
    }

    pub fn run(self: *Self) !void {
        while (!self.exit) {
            try switch (self.params.words) {
                inline 1...bignum.MAXWORDS => |v| self.renderForParams(bignum.BigInt(64 * v)),
                else => self.renderForParams(bignum.BigInt(64 * bignum.MAXWORDS)),
            };

            if (self.pause) {
                while (!self.exit) {
                    if (try self.handleInput())
                        break;
                }
            }
        }
    }

    pub fn deinit(self: *Self) void {
        self.storage.deinit();
        sdl2.SDL_DestroyRenderer(self.renderer);
        sdl2.SDL_DestroyWindow(self.window);
    }

    fn renderBlock(self: Self, bpx: u32, bpy: u32, block: *mandel.BlockData) void {
        var rect : sdl2.SDL_FRect = undefined;

        const rectSize: u32 = @as(u32, 1) << self.params.magShift;
        var oy: u32 = 0;
        while (oy < block.sz) : (oy += rectSize) {
            const py = bpy + oy;
            var ox: u32 = 0;
            while (ox < block.sz) : (ox += rectSize) {
                const px = bpx + ox;

                const current = block.iter(ox, oy);

                var g: u8 = @intCast(u8, if (current < 0 or current > self.params.iters) 0 else current & 0xff);
                // const a = @intCast(u8, if (current <= 0) 0x0 else if (current > self.params.iters) 0xff else @intCast(u32, (@intCast(u32, current) * 255 / self.params.iters)) & 0xff);
                const a : u8 = 0xff;
                g = (g << 7) | (g >> 1);
                _ = sdl2.SDL_SetRenderDrawColor(self.renderer, g, g, g, a);

                rect.x = @intToFloat(f32, px);
                rect.y = @intToFloat(f32, py);
                rect.w = @intToFloat(f32, rectSize);
                rect.h = @intToFloat(f32, rectSize);
                _ = sdl2.SDL_RenderFillRectF(self.renderer, &rect);
            }
        }
    }

    fn clearParamsT(self: *Self, calcSize : u32, comptime BigIntType : type) !void {
        var workLoad = try Params.BlockWorkMaker(BigIntType).init(self.alloc, self.params, &self.storage,
            calcSize, 0, 0,
        );
        defer workLoad.deinit(self.alloc);

        for (workLoad.blocks) |block| {
            std.mem.set(i32, block.data.iters, 0);
        }
    }

    fn clearParams(self: *Self, calcSize : u32) !void {
        try switch (self.params.words) {
            inline 1...bignum.MAXWORDS => |v| self.clearParamsT(calcSize, bignum.BigInt(64 * v)),
            else => self.clearParamsT(calcSize, bignum.BigInt(64 * bignum.MAXWORDS)),
        };
    }

    fn renderForParams(self: *Self, comptime BigIntType: type) !void {
        std.debug.print("Rendering with BigInt({s}) at x= {s}, y= {s}, zoom= {}, iters= {}...\n", .{ @typeName(BigIntType), self.params.cx, self.params.cy, self.params.zoom, self.params.iters });

        _ = sdl2.SDL_SetRenderDrawColor(self.renderer, 0, 0, 0, 0);
        _ = sdl2.SDL_RenderClear(self.renderer);

        var timer = try std.time.Timer.start();
        const timeStart = timer.lap();

        var wx : c_int = undefined;
        var wy : c_int = undefined;
        _ = sdl2.SDL_GetWindowSize(self.window, &wx, &wy);
        const calcSize = @intCast(u32, std.math.max(wx, wy));

        var workLoad = try Params.BlockWorkMaker(BigIntType).init(self.alloc, self.params, &self.storage,
                calcSize, 0, 0,
        );
        defer workLoad.deinit(self.alloc);

        std.debug.print("time for making blocks = {}\n", .{timer.lap() - timeStart});

        for (workLoad.blocks) |block| {
            self.renderBlock(block.px, block.py, block.data);
        }

        _ = sdl2.SDL_RenderPresent(self.renderer);

        std.debug.print("time for first block draw = {}\n", .{timer.lap() - timeStart});

        if (self.pause) {
            return;
        }

        const BlockType = @TypeOf(workLoad.blocks[0]);

        const ThreadInfo = struct {
            fn calcThread(stop: *std.atomic.Atomic(bool), readyBlocks: *std.ArrayList(*BlockType), readyMutex: *std.Thread.Mutex,
                    blockIndex: *std.atomic.Atomic(usize), blocks: []BlockType,
                    blocksLeft: *std.atomic.Atomic(usize)) !void {
                while (!stop.load(.Unordered)) {
                    var myIndex = blockIndex.fetchAdd(1, .SeqCst);
                    if (myIndex >= blocks.len) break;

                    var block = &blocks[myIndex];

                    try block.calculate(stop);

                    readyMutex.lock();
                    defer readyMutex.unlock();
                    try readyBlocks.append(block);

                    _ = blocksLeft.fetchSub(1, .SeqCst);
                }
            }
        };

        // make worker threads
        var stop = std.atomic.Atomic(bool).init(false);
        var readyBlocks = std.ArrayList(*BlockType).init(self.alloc);
        defer readyBlocks.deinit();
        var readyMutex: std.Thread.Mutex = .{};
        var blockIndex = std.atomic.Atomic(usize).init(0);
        var blocksLeft = std.atomic.Atomic(usize).init(workLoad.blocks.len);

        const NTHREADS: usize = std.math.min(workLoad.blocks.len, std.math.min(config.MAXTHREADS, try std.Thread.getCpuCount()));
        var threads: [256]std.Thread = undefined;

        var nt: u32 = 0;
        while (nt < NTHREADS) : (nt += 1) {
            threads[nt] = try std.Thread.spawn(.{}, ThreadInfo.calcThread, .{
                &stop, &readyBlocks, &readyMutex,
                &blockIndex, workLoad.blocks, &blocksLeft
            });
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

            var block: ?*BlockType = null;
            var any = false;
            {
                readyMutex.lock();
                defer readyMutex.unlock();
                while (true) {
                    block = readyBlocks.popOrNull();
                    if (block == null) {
                        if (blocksLeft.load(.Acquire) == 0) {
                            done = true;
                            std.debug.print("time for calculation   = {d}s\n", .{@intToFloat(f64, (timer.lap() - timeStart)) / 1.0e9});
                        }
                        break;
                    }
                    self.renderBlock(block.?.px, block.?.py, block.?.data);
                    any = true;
                }
            }

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

    fn setPause(self: *Self, p: bool) void {
        std.debug.print("{s}\n", .{if (p) "PAUSE" else "RESUME"});
        self.pause = p;
    }

    fn handleInput(self: *Self) !bool {
        var wx : c_int = undefined;
        var wy : c_int = undefined;
        _ = sdl2.SDL_GetWindowSize(self.window, &wx, &wy);
        const calcSize = @intCast(u32, std.math.max(wx, wy));
        const blockSize = self.storage.blockSize;
        var ps = &self.params;
        var event: sdl2.SDL_Event = undefined;
        var cont = false;
        while (sdl2.SDL_WaitEventTimeout(&event, 10) != 0) {
            var words = ps.words;
            var span = ps.span(calcSize);
            if (event.type == sdl2.SDL_QUIT) {
                self.exit = true;
                std.debug.print("Cancelling...\n", .{});
                cont = true;
            } else if (event.type == sdl2.SDL_KEYDOWN) {
                var ke = @ptrCast(*sdl2.SDL_KeyboardEvent, &event);
                cont = true;

                var shift = (ke.state & sdl2.KMOD_SHIFT) != 0;
                var minShift : u5 = std.math.min(1, @intCast(u5, std.math.log2_int(u32, @intCast(u32, std.math.min(wx, wy)) / blockSize)));
                const shiftAmt: u5 = if ((ke.keysym.mod & sdl2.KMOD_LSHIFT) != 0) (minShift + 3) else minShift;
                std.debug.print("Min shift= {}, shift = {}, span = {}\n", .{ minShift, shiftAmt, @floatCast(f64, span) });
                _ = switch (ke.keysym.sym) {
                    sdl2.SDLK_ESCAPE => { self.exit = true; std.debug.print("Cancelling...\n", .{}); },
                    sdl2.SDLK_PLUS, sdl2.SDLK_EQUALS => ps.zoom += 1,
                    sdl2.SDLK_MINUS => ps.zoom = if (ps.zoom > 1) ps.zoom - 1 else 0,
                    sdl2.SDLK_UP => try bignum.faddShift(self.alloc, &ps.cy, words, -span, shiftAmt),
                    sdl2.SDLK_DOWN => try bignum.faddShift(self.alloc, &ps.cy, words, span, shiftAmt),
                    sdl2.SDLK_LEFT => try bignum.faddShift(self.alloc, &ps.cx, words, -span, shiftAmt),
                    sdl2.SDLK_RIGHT => try bignum.faddShift(self.alloc, &ps.cx, words, span, shiftAmt),
                    sdl2.SDLK_PAGEDOWN => ps.iters = @intCast(u32, std.math.max(0, @intCast(i32, ps.iters) - if (shift) @intCast(i32, 250) else @intCast(i32, 50))),
                    sdl2.SDLK_PAGEUP => ps.iters += if (shift) @intCast(u32, 250) else @intCast(u32, 50),
                    sdl2.SDLK_LEFTBRACKET => ps.words = if (words > 1) words - 1 else 1,
                    sdl2.SDLK_RIGHTBRACKET => ps.words += 1,
                    sdl2.SDLK_SPACE => {
                        try bignum.falign(self.alloc, &ps.cx, words);
                        try bignum.falign(self.alloc, &ps.cy, words);
                    },
                    sdl2.SDLK_r => try self.clearParams(calcSize),
                    sdl2.SDLK_F5 => {
                        try ps.writeConfig(self.alloc, self.saveFile);
                        cont = false;
                    },
                    sdl2.SDLK_F9 => {
                        if (ps.readConfig(self.alloc, self.saveFile)) {} else |err| {
                            std.debug.print("Failed to load: {}\n", .{err});
                            cont = false;
                        }
                    },
                    sdl2.SDLK_F4 => {
                        var file = try files.RenderedFile.init(&self.params, &self.storage);
                        // if (file.save(try std.fmt.allocPrint(self.alloc, "{s}.tst", .{self.dataFile}))) {
                        if (file.save(self.dataFile)) {
                            cont = true;
                        } else |err| {
                            std.debug.print("Failed to save: {}\n", .{err});
                            cont = false;
                        }
                    },
                    sdl2.SDLK_F8 => {
                        var file = try files.RenderedFile.init(&self.params, &self.storage);
                        if (file.load(self.dataFile)) {
                            cont = true;
                        } else |err| {
                            std.debug.print("Failed to load: {}\n", .{err});
                            cont = false;
                        }
                    },
                    sdl2.SDLK_PAUSE => { self.setPause(!self.pause); cont = true; },
                    '0'...'9' => |v| ps.magShift = @intCast(u5, v - '0'),
                    else => cont = false,
                };
            } else if (event.type == sdl2.SDL_MOUSEBUTTONDOWN) {
                var me = @ptrCast(*sdl2.SDL_MouseButtonEvent, &event);
                const shift = (sdl2.SDL_GetModState() & sdl2.KMOD_SHIFT) != 0;
                if (me.button == 1 and shift) {
                    cont = true;
                    // recenter on mouse position
                    try bignum.fadj(self.alloc, &ps.cy, words, span, (me.y + (blockSize >> 1)) & ~(blockSize - 1), @intCast(i32, ps.sy));
                    try bignum.fadj(self.alloc, &ps.cx, words, span, (me.x + (blockSize >> 1)) & ~(blockSize - 1), @intCast(i32, ps.sx));
                }
            }
        }
        return cont;
    }
};
