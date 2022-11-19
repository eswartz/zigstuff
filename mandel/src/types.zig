
pub const XY = struct {
    x: u32, y: u32,

    pub fn hash(a: @This()) u64 {
        return a.x ^ a.y;
    }
    pub fn eql(me: @This(), other: @This()) bool {
        return me.x == other.x and me.y == other.y;
    }
};
