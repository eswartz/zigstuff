

pub const sdl2 = @cImport({
    @cInclude("SDL2/SDL.h");
});

pub const gmp = @cImport({
    @cInclude("gmp.h");
});

pub const zlib = @cImport({
    @cInclude("zlib.h");
});
