
set(PROGRAMS optimisers)

foreach(prog ${PROGRAMS})
    add_executable(${prog} ${prog}.cpp)
    target_link_libraries(${prog} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ook)
endforeach(prog)
