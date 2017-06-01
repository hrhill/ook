set(PROGRAMS version
             nlgc_comparison)

foreach(ex ${PROGRAMS})
    add_executable(${ex} ${ex}.cpp)
    target_link_libraries(${ex} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ook)
endforeach(ex)


add_subdirectory(validated_numerics)
add_subdirectory(optimisers)

