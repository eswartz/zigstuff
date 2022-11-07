set history save on
set breakpoint pending on
catch throw
set height 0
set confirm off

#set scheduler-locking step
set print elements 100
set print thread-events off
set print inferior-events off

set print asm-demangle on
set print object on
set print static-members off

## used by threading
handle SIGPWR nostop noprint
handle SIGXCPU nostop noprint

define rr
    dir 
    r
end

define jt
    tb $arg0
    jump $arg0
end

