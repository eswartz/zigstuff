const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;

const config = @import("config.zig");

const bignum = @import("bignum.zig");
const Mandel = bignum.Mandel;

const XY = @import("types.zig").XY;

const algo = @import("algo.zig");
const MandelParams = algo.Params;
const cache = @import("cache.zig");
const MandelStorage = cache.MandelStorage;

const cimports = @import("cimports.zig");
const sdl2 = cimports.sdl2;

const files = @import("files.zig");

pub const Viewer = struct {
    const Self = @This();

    alloc: Allocator,
    window: *sdl2.SDL_Window,
    renderer: *sdl2.SDL_Renderer,
    params: MandelParams,

    currentBlockList: *MandelParams.BlockList,
    exit: bool,
    pause: bool,
    saveFile : []u8,
    dataFile : []u8,
    storage : MandelStorage,
    updateMetadata : bool,
    autoSave : bool,
    autoFrom : u16,
    autoTo : u16,

    pub fn init(alloc: Allocator) !Self {
        var winSize : c_int = 1024;
        var frameSize : u32 = 1024;
        var blockSize : u16 = 128;

        var updateMetadata : bool = false;
        var autoFrom : u16 = 0;
        var autoTo : u16 = 0;
        var saveName : []u8 = "";
        var dataName : []u8 = "";
        {
            var argIter = std.process.args();
            _ = argIter.next(); // name
            while (argIter.next()) |arg| {
                if (std.mem.eql(u8, arg, "-h")) {
                    std.debug.print(
                        \\ help: [options] .json ...%Z....dat
                        \\
                        \\ -s pot: render size
                        \\ -r pot: window size
                        \\ -b pot: block size
                        \\ -a N..M: automatically render zoom levels N..M
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
                } else if (std.mem.eql(u8, arg, "-a")) {
                    var argVal = argIter.next().?;
                    if (std.mem.indexOf(u8, argVal, "..")) |idx| {
                        autoFrom = try std.fmt.parseInt(u16, argVal[0..idx], 10);
                        autoTo = try std.fmt.parseInt(u16, argVal[idx+2..argVal.len], 10);
                        if (autoFrom > autoTo) return error.RangeBackwards;
                    } else {
                        autoFrom = try std.fmt.parseInt(u16, argVal, 10);
                        autoTo = autoFrom + 1;
                    }
                } else if (std.mem.eql(u8, arg, "-u")) {
                    updateMetadata = true;
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

        if (updateMetadata and autoFrom >= autoTo) {
            return error.NeedAutosave;
        }

        if (blockSize > winSize) blockSize = @intCast(u16, winSize);
        if (winSize > frameSize) frameSize = @intCast(u32, winSize);

        if (saveName.len == 0) saveName = try alloc.dupe(u8, "save.json");
        if (dataName.len == 0) dataName = try alloc.dupe(u8, "data/zoom10_000.dat");

        var store = try MandelStorage.init(alloc, blockSize);

        var params = try MandelParams.init(alloc);

        params.sx = frameSize;
        params.sy = frameSize;
        if (autoFrom < autoTo) {
            params.zoom = autoFrom;
            params.magShift = 0;
        }

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
            .currentBlockList = undefined,
            .exit = false,
            .pause = false,
            .saveFile = saveName,
            .dataFile = dataName,
            .storage = store,
            .updateMetadata = updateMetadata,
            .autoSave = autoFrom < autoTo,
            .autoFrom = autoFrom,
            .autoTo = autoTo,
        };
    }

    pub fn run(self: *Self) !void {
        if (self.autoSave) {
            try files.JsonFormat.loadConfig(&self.params, self.alloc, self.saveFile);
            var rootParams = self.params;
            rootParams.cx = try self.alloc.dupe(u8, self.params.cx);
            rootParams.cy = try self.alloc.dupe(u8, self.params.cy);
            self.params.zoom = self.autoFrom;
            self.params.magShift = 0;
            self.params.words = self.params.getIntSizeForZoom();

            while (true) {
                // load any previous data
                if (try self.loadData()) |_| {
                    // restore guiding coords
                    self.params.setCx(self.alloc, try self.alloc.dupe(u8, rootParams.cx));
                    self.params.setCy(self.alloc, try self.alloc.dupe(u8, rootParams.cy));
                    try self.render();
                }

                const changed = try self.recalc();
                std.debug.print("changed={}\n", .{changed});

                while (true) {
                    const cont = try self.handleInput();
                    if (self.exit)
                        return;
                    if (!self.pause)
                        break;
                    if (cont) {
                        try self.render();
                    }
                }

                if (self.pause) {
                    std.debug.print("Paused, not saving\n", .{});
                    continue;
                } if (!self.updateMetadata and !changed) {
                    std.debug.print("No changes for {} of {}\n", .{ self.params.zoom, self.autoTo });
                } else {
                    std.debug.print("Saving {} of {}\n", .{ self.params.zoom, self.autoTo });
                    if (self.saveData()) |_| {
                        // save memory
                    } else |err| {
                        return err;
                    }
                }

                self.params.setZoom(self.params.zoom + 1);
                if (self.params.zoom > self.autoTo) {
                    self.exit = true;
                    return;
                }

            } else |err| {
                return err;
            }
        } else {
            if (files.JsonFormat.loadConfig(&self.params, self.alloc, self.saveFile)) {} else |err| {
                std.debug.print("Failed to load: {}\n", .{err});
            }

            while (!self.exit) {
                if (!self.pause) {
                    if (try self.recalc()) {
                        continue;
                    }
                } else {
                    try self.render();
                }

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

    fn renderBlock(self: Self, bpx: u32, bpy: u32, block: *cache.BlockData) void {
        const rectSize: u32 = @as(u32, 1) << self.params.magShift;

        // predetermine if it's complete
        var isComplete = true;
        {
            var oy: u32 = 0;
            loop: while (oy < block.sz) : (oy += rectSize) {
                var ox: u32 = 0;
                while (ox < block.sz) : (ox += rectSize) {
                    const current = block.iter(ox, oy);
                    if (current <= 0 and -current < self.params.iters) {
                        isComplete = false;
                        break :loop;
                    }
                }
            }
        }

        var rect : sdl2.SDL_Rect = .{ .x = @intCast(c_int, bpx), .y = @intCast(c_int, bpy), .w = block.sz, .h = block.sz };
        {
            // clear on-screen block
            _ = sdl2.SDL_SetRenderDrawColor(self.renderer, 0, 0, 0, 0);
            _ = sdl2.SDL_RenderFillRect(self.renderer, &rect);
        }

        const alpha : u8 = if (isComplete) 0xff else 0x40;
        if (!isComplete) {
            _ = sdl2.SDL_SetRenderDrawBlendMode(self.renderer, sdl2.SDL_BLENDMODE_BLEND);
        } else {
            _ = sdl2.SDL_SetRenderDrawBlendMode(self.renderer, sdl2.SDL_BLENDMODE_NONE);
        }

        // generate pixels as rects (per magShift)
        rect.w = @intCast(c_int, rectSize);
        rect.h = @intCast(c_int, rectSize);

        var oy: u32 = 0;
        while (oy < block.sz) : (oy += rectSize) {
            const py = bpy + oy;
            var ox: u32 = 0;
            while (ox < block.sz) : (ox += rectSize) {
                const px = bpx + ox;

                const current = block.iter(ox, oy);

                var gamma: u8 = @intCast(u8, if (current < 0 or current > self.params.iters) 0 else current & 0xff);
                // const a = @intCast(u8, if (current <= 0) 0x0 else if (current > self.params.iters) 0xff else @intCast(u32, (@intCast(u32, current) * 255 / self.params.iters)) & 0xff);
                const a : u8 = alpha;
                const r = (gamma << 5) | (gamma >> 3);
                const g = (gamma << 6) | (gamma >> 2);
                const b = (gamma << 7) | (gamma >> 1);
                _ = sdl2.SDL_SetRenderDrawColor(self.renderer, r, g, b, a);

                rect.x = @intCast(c_int, px);
                rect.y = @intCast(c_int, py);
                _ = sdl2.SDL_RenderFillRect(self.renderer, &rect);
            }
        }

        _ = sdl2.SDL_SetRenderDrawBlendMode(self.renderer, sdl2.SDL_BLENDMODE_NONE);
    }

    fn clearParamsT(self: *Self, comptime BigIntType : type, zoomOffs : i16) !void {
        const calcSize = self.getCalcSize();
        var params = self.params;
        if (params.zoom == 0 and zoomOffs < 0) return;

        params.zoom = @intCast(u16, @intCast(i16, params.zoom) + zoomOffs);
        std.debug.print("Resetting blocks at zoom {}\n", .{ params.zoom });
        var blockList = try MandelParams.BlockGenerator(BigIntType).init(self.alloc, params, &self.storage, calcSize);
        defer blockList.deinit(self.alloc);

        for (blockList.blocks) |block| {
            std.mem.set(i32, block.data.iters, 0);
        }
    }

    fn clearParams(self: *Self) !void {
        var i : i16 = -1;
        while (i < 2) : (i += 1) {
            // try switch (self.params.words) {
            //     inline 1...config.MAXWORDS => |v| self.clearParamsT(bignum.BigInt(64 * v), i),
            //     else => self.clearParamsT(bignum.BigInt(64 * config.MAXWORDS), i),
            // };
            try self.clearParamsT(bignum.BigInt(64 * config.MAXWORDS), i);
        }
    }

    fn getCalcSize(self: Self) u32 {
        var wx : c_int = undefined;
        var wy : c_int = undefined;
        _ = sdl2.SDL_GetWindowSize(self.window, &wx, &wy);
        const calcSize = @intCast(u32, std.math.floorPowerOfTwo(c_int, std.math.max(wx, wy)));
        return calcSize;
    }

    // fn ensureBlockList(self: *Self, comptime BigIntType: type) !BlockList {
    //     if (self.currentBlockList)
    //     var blockList = try MandelParams.BlockGenerator(BigIntType).init(self.alloc, self.params, &self.storage, self.getCalcSize());
    //     defer blockList.deinit(self.alloc);

    // }

    fn renderT(self: *Self, comptime BigIntType: type) !void {
        var blockList = try MandelParams.BlockGenerator(BigIntType).init(self.alloc, self.params, &self.storage, self.getCalcSize());
        defer blockList.deinit(self.alloc);

        _ = sdl2.SDL_SetRenderDrawColor(self.renderer, 0, 0, 0, 0);
        _ = sdl2.SDL_RenderClear(self.renderer);

        for (blockList.blocks) |block| {
            self.renderBlock(block.px, block.py, block.data);
        }

        _ = sdl2.SDL_SetRenderDrawColor(self.renderer, 255, 255, 255, 32);

        const calcSizeH = @intCast(c_int, self.getCalcSize()) >> 1;
        _ = sdl2.SDL_RenderDrawLine(self.renderer, calcSizeH-8, calcSizeH, calcSizeH+8, calcSizeH);
        _ = sdl2.SDL_RenderDrawLine(self.renderer, calcSizeH, calcSizeH-8, calcSizeH, calcSizeH+8);

        _ = sdl2.SDL_RenderPresent(self.renderer);
    }

    fn render(self: *Self) !void {
        try switch (self.params.words) {
            inline 1...bignum.MAXWORDS => |v| self.renderT(bignum.BigInt(64 * v)),
            else => self.renderT(bignum.BigInt(64 * config.MAXWORDS)),
        };
    }

    fn recalcT(self: *Self, comptime BigIntType: type) !bool {
        std.debug.print("Rendering with BigInt({s}) at x= {s}, y= {s}, zoom= {}, iters= {}...\n", .{ @typeName(BigIntType), self.params.cx, self.params.cy, self.params.zoom, self.params.iters });

        var timer = try std.time.Timer.start();
        const timeStart = timer.lap();

        const calcSize = self.getCalcSize();

        var blockList = try MandelParams.BlockGenerator(BigIntType).init(self.alloc, self.params, &self.storage, calcSize);
        defer blockList.deinit(self.alloc);

        std.debug.print("time for making blocks = {}\n", .{timer.lap() - timeStart});

        const BlockType = @TypeOf(blockList.blocks[0]);

        const ThreadInfo = struct {
            fn calcThread(stop: *std.atomic.Atomic(bool),
                    readyBlocks: *std.ArrayList(*BlockType), readyMutex: *std.Thread.Mutex,
                    blockIndex: *std.atomic.Atomic(usize), blocks: []BlockType,
                    blocksLeft: *std.atomic.Atomic(usize),
                    blocksChanged: *std.atomic.Atomic(usize),
                    ) !void {
                while (!stop.load(.Unordered)) {
                    var myIndex = blockIndex.fetchAdd(1, .SeqCst);
                    if (myIndex >= blocks.len) break;

                    var block = &blocks[myIndex];

                    if (try block.calculate(stop)) {
                        _ = blocksChanged.fetchAdd(1, .SeqCst);
                    }

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
        var blocksChanged = std.atomic.Atomic(usize).init(0);
        var blockIndex = std.atomic.Atomic(usize).init(0);
        var blocksLeft = std.atomic.Atomic(usize).init(blockList.blocks.len);

        const NTHREADS: usize = 3 * std.math.min(blockList.blocks.len, std.math.min(config.MAXTHREADS, try std.Thread.getCpuCount())) / 4;
        var threads: [256]std.Thread = undefined;

        var nt: u32 = 0;
        while (nt < NTHREADS) : (nt += 1) {
            threads[nt] = try std.Thread.spawn(.{}, ThreadInfo.calcThread, .{
                &stop, &readyBlocks, &readyMutex,
                &blockIndex, blockList.blocks, &blocksLeft, &blocksChanged,
            });
        }

        std.debug.print("time for thread launch = {}\n", .{timer.lap() - timeStart});

        // render as work arrives
        var done: bool = false;

        while (!done) {
            if (try self.handleInput()) {
                // params changed, abort recalculation
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
                        if (blocksLeft.load(.SeqCst) == 0) {
                            done = true;
                            std.debug.print("time for calculation   = {d}s\n", .{@intToFloat(f64, (timer.lap() - timeStart)) / 1.0e9});
                        }
                        break;
                    }
                    self.renderBlock(block.?.px, block.?.py, block.?.data);
                    any = true;
                }
            }

            // render anything that updated
            if (any) {
                _ = sdl2.SDL_RenderPresent(self.renderer);
            }
        }

        // clean up
        stop.store(true, .Unordered);
        for (threads[0..NTHREADS]) |thread| {
            thread.join();
        }

        return blocksChanged.load(.SeqCst) > 0;
    }

    fn recalc(self: *Self) !bool {
        return try switch (self.params.words) {
            inline 1...config.MAXWORDS => |v| self.recalcT(bignum.BigInt(64 * v)),
            else => self.recalcT(bignum.BigInt(64 * config.MAXWORDS)),
        };
    }

    fn setPause(self: *Self, p: bool) void {
        std.debug.print("{s}\n", .{if (p) "PAUSE" else "RESUME"});
        self.pause = p;
    }

    fn formatDataFile(self: Self) ![]u8 {
        var zoomStr = try std.fmt.allocPrint(self.alloc, "{:3}", .{ self.params.zoom }); defer self.alloc.free(zoomStr);
        const replSize = std.mem.replacementSize(u8, self.dataFile, "%Z", zoomStr);
        const data = try self.alloc.alloc(u8, replSize);
        _ = std.mem.replace(u8, self.dataFile, "%Z", zoomStr, data);
        _ = std.mem.replace(u8, data, " ", "0", data);
        std.debug.print("Expanded {s} to {s}\n", .{self.dataFile, data});
        return data;
    }

    fn saveData(self: *Self) !bool {
        var data = try self.formatDataFile(); defer self.alloc.free(data);
        var file = try files.RenderedFile.init(&self.params, &self.storage);
        if (file.save(data)) {
            return true;
        } else |err| {
            std.debug.print("Failed to save: {}\n", .{err});
            return false;
        }
    }

    fn loadData(self: *Self) !bool {
        var data = try self.formatDataFile(); defer self.alloc.free(data);
        var file = try files.RenderedFile.init(&self.params, &self.storage);
        if (file.load(data)) |complete| {
            return complete;
        } else |err| {
            std.debug.print("Failed to load: {}\n", .{err});
            return false;
        }
    }

    /// Look for an interesting area nearby where there's a hole, defined
    /// as an area whose Block has the most number of high-value unfinished iterations
    fn findHoleT(self: *Self, comptime BigIntType: type) !XY {
        var blockList = try MandelParams.BlockGenerator(BigIntType).init(self.storage.alloc, self.params, &self.storage, self.params.sx);
        defer blockList.deinit(self.storage.alloc);

        const ps = self.storage.blockSize;

        var deepestIter : u32 = 0;
        var deepestNumUnset : u32 = 0;
        var deepestScore : f64 = 0;
        var deepestPx : u32 = 0;
        var deepestPy : u32 = 0;

        const mag: u32 = @as(u32, 1) << self.params.magShift;
        for (blockList.blocks) |block| {
            var numUnset : u32 = 0;
            var maxIter : u32 = 0;
            var score : f64 = 0;

            var oy : u32 = 0;
            while (oy < ps) : (oy += mag) {
                var ox : u32 = 0;
                while (ox < ps) : (ox += mag) {
                    const iter = block.data.iter(ox, oy);
                    if (iter < 0 or iter >= self.params.iters) numUnset += 1;
                    const absIter = @intCast(u32, if (iter < 0) -iter else iter);
                    score += @intToFloat(f64, absIter);
                    if (absIter > maxIter) {
                        maxIter = absIter;
                    }
                }
            }
            if (score > deepestScore and (numUnset > deepestNumUnset or maxIter >= deepestIter)) {
                deepestPx = block.px;
                deepestPy = block.py;
                deepestIter = maxIter;
                deepestScore = score;
                std.debug.print("--> {},{} = {} #{}\n", .{ block.px, block.py, maxIter, numUnset });
            }
        }
        if (deepestIter == 0) return error.NoDeepestFound;
        return XY{ .x = deepestPx, .y = deepestPy };
    }

    fn findHole(self: *Self) !XY {
        return try switch (self.params.words) {
            inline 1...config.MAXWORDS => |v| self.findHoleT(bignum.BigInt(64 * v)),
            else => self.findHoleT(bignum.BigInt(64 * config.MAXWORDS)),
        };
    }

    fn handleInput(self: *Self) !bool {
        const blockSize = self.storage.blockSize;

        var minShift : u5 = std.math.min(1, @intCast(u5, std.math.log2_int(u32, @intCast(u32, self.getCalcSize()) / blockSize)));

        var ps = &self.params;
        var event: sdl2.SDL_Event = undefined;
        var cont = false;

        while (sdl2.SDL_WaitEventTimeout(&event, 10) != 0) {
            var words = ps.words;
            var span = ps.span(self.getCalcSize());
            switch (event.type) {
            sdl2.SDL_QUIT => {
                self.exit = true;
                std.debug.print("Cancelling...\n", .{});
                cont = true;
            },
            sdl2.SDL_WINDOWEVENT => {
                var we = @ptrCast(*sdl2.SDL_WindowEvent, &event);
                switch (we.event) {
                sdl2.SDL_WINDOWEVENT_EXPOSED, sdl2.SDL_WINDOWEVENT_RESIZED, sdl2.SDL_WINDOWEVENT_SHOWN => |e| {
                    std.debug.print("Redrawing {}\n", .{e});
                    try self.render();
                },
                else => {}
                }
            },
            sdl2.SDL_KEYDOWN => {
                var ke = @ptrCast(*sdl2.SDL_KeyboardEvent, &event);
                const shiftAmt: u5 = if ((ke.keysym.mod & sdl2.KMOD_LSHIFT) != 0) (minShift + 3) else minShift;

                // std.debug.print("Min shift= {}, shift = {}, span = {}\n", .{ minShift, shiftAmt, @floatCast(f64, span) });
                cont = !self.autoSave;

                var shift = (ke.keysym.mod & sdl2.KMOD_SHIFT) != 0;
                var ctrl = (ke.keysym.mod & sdl2.KMOD_CTRL) != 0;

                _ = switch (ke.keysym.sym) {
                    sdl2.SDLK_ESCAPE => { self.exit = true; std.debug.print("Cancelling...\n", .{}); },
                    sdl2.SDLK_F1, sdl2.SDLK_SLASH, sdl2.SDLK_QUESTION => {
                        std.debug.print(
                            \\ Help:
                            \\   ?/F1: help
                            \\   Esc: exit
                            \\   Pause: pause rendering and autosaving
                            \\   Space: re-render without calculation
                            \\   +: zoom in
                            \\   -: zoom out
                            \\   pgup: +50 iterations (shift: +250)
                            \\   pgdn: -50 iterations (shift: -250)
                            \\   up/down/left/right: scroll view half a page, shift+: scroll 1/8 page
                            \\   0-9: sample 2^n pixels
                            \\   [: -64 bits precision
                            \\   ]: +64 bits precision
                            \\   \: save zoom -> precision map
                            \\   ctrl-\: save zoom -> iters map
                            \\   ctrl-a: align view to precision
                            \\   ctrl-r: clear render at (-1, 0, 1) zoom levels
                            \\   F5: save .json params
                            \\   F9: load .json params
                            \\   F4: save .dat file
                            \\   F8: load .dat file
                            \\   F12: search for a hole and move there
                            \\
                            , .{});
                        cont = false;
                    },
                    sdl2.SDLK_SPACE => try self.render(),
                    sdl2.SDLK_RETURN => cont = true,
                    sdl2.SDLK_PLUS, sdl2.SDLK_EQUALS => ps.setZoomInteractive(ps.zoom + 1),
                    sdl2.SDLK_MINUS => ps.setZoomInteractive(if (ps.zoom > 1) ps.zoom - 1 else 0),
                    sdl2.SDLK_UP => if (self.pause or !self.autoSave) try bignum.faddShift(self.alloc, &ps.cy, words, -span, shiftAmt),
                    sdl2.SDLK_DOWN => if (self.pause or !self.autoSave) try bignum.faddShift(self.alloc, &ps.cy, words, span, shiftAmt),
                    sdl2.SDLK_LEFT => if (self.pause or !self.autoSave) try bignum.faddShift(self.alloc, &ps.cx, words, -span, shiftAmt),
                    sdl2.SDLK_RIGHT => if (self.pause or !self.autoSave) try bignum.faddShift(self.alloc, &ps.cx, words, span, shiftAmt),
                    sdl2.SDLK_PAGEDOWN => ps.iters = @intCast(u32, std.math.max(0, @intCast(i32, ps.iters) - (if (shift) @intCast(i32, 250) else @intCast(i32, 50)))),
                    sdl2.SDLK_PAGEUP => ps.iters += if (shift) @intCast(u32, 250) else @intCast(u32, 50),
                    sdl2.SDLK_LEFTBRACKET => ps.words = if (words > 1) words - 1 else 1,
                    sdl2.SDLK_RIGHTBRACKET => ps.words = if (words < bignum.MAXWORDS) words + 1 else bignum.MAXWORDS,
                    sdl2.SDLK_BACKSLASH => { if (ctrl) try ps.recordZoomIters() else try ps.recordZoomBits(); cont = false; },
                    sdl2.SDLK_a => if (ctrl) {
                        try bignum.falign(self.alloc, &ps.cx, words);
                        try bignum.falign(self.alloc, &ps.cy, words);
                    },
                    sdl2.SDLK_r => if (ctrl) try self.clearParams(),
                    sdl2.SDLK_F5 => {
                        try files.JsonFormat.writeConfig(ps.*, self.alloc, self.saveFile);
                        cont = false;
                    },
                    sdl2.SDLK_F9 => {
                        if (files.JsonFormat.loadConfig(ps, self.alloc, self.saveFile)) {} else |err| {
                            std.debug.print("Failed to load: {}\n", .{err});
                            cont = false;
                        }
                    },
                    sdl2.SDLK_F4 => cont = try self.saveData(),
                    sdl2.SDLK_F8 => { cont = try self.loadData(); try self.render(); },
                    sdl2.SDLK_PAUSE => { self.setPause(!self.pause); cont = true; },
                    '0'...'9' => |v| ps.magShift = @intCast(u5, v - '0'),
                    sdl2.SDLK_F12 => {
                        if (self.findHole()) |xy| {
                            std.debug.print("xy = {}, {}\n", .{xy.x, xy.y});
                            try bignum.fadj(self.alloc, &ps.cx, words, span, @intCast(i32, xy.x), @intCast(i32, ps.sx));
                            try bignum.fadj(self.alloc, &ps.cy, words, span, @intCast(i32, xy.y), @intCast(i32, ps.sy));
                            // cont  = false;
                        } else |err| if (err == error.NoDeepestFound) {
                            // fine
                            cont  = false;
                        } else {
                            return err;
                        }
                    },
                    else => cont = false,
                };
            },
            sdl2.SDL_MOUSEBUTTONDOWN => {
                var me = @ptrCast(*sdl2.SDL_MouseButtonEvent, &event);
                const shift = (sdl2.SDL_GetModState() & sdl2.KMOD_SHIFT) != 0;
                if (me.button == 1 and shift) {
                    cont = true;
                    // recenter on mouse position
                    try bignum.fadj(self.alloc, &ps.cx, words, span, (me.x + (blockSize >> 1)) & ~(blockSize - 1), @intCast(i32, ps.sx));
                    try bignum.fadj(self.alloc, &ps.cy, words, span, (me.y + (blockSize >> 1)) & ~(blockSize - 1), @intCast(i32, ps.sy));
                }
            },
            else => {}
            }
        }
        return cont;
    }
};
