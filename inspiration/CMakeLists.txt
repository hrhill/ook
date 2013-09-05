set(PROGRAMS line_search_test)

foreach(prog ${PROGRAMS})
    add_executable(${prog} ${prog}.cpp dcsrch.cpp dcstep.cpp)
endforeach(prog)
