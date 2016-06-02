
set(PROGRAMS optimisers)

foreach(prog ${PROGRAMS})
    add_executable(${prog} ${prog}.cpp)
    target_link_libraries(${prog} ${Boost_PROGRAM_OPTIONS_LIBRARY} ook)
endforeach(prog)
