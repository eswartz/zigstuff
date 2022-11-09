
pub const MAXWORDS = 16;
pub const MAXBITS = MAXWORDS << 6;
pub const BigIntMax = @Type(.{ .Int = .{ .signedness = .signed, .bits = MAXBITS } });
pub const MAXTHREADS = 256;
