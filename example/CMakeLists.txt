set(PROGRAMS version
             nlgc_comparison)

foreach(ex ${PROGRAMS})
    add_executable(${ex} ${ex}.cpp)
    target_link_libraries(${ex} ${Boost_PROGRAM_OPTIONS_LIBRARY} ook)
endforeach(ex)


add_subdirectory(validated_numerics)
add_subdirectory(optimisers)

