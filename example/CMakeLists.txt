set(PROGRAMS version)

foreach(ex ${PROGRAMS})
    add_executable(${ex} ${ex}.cpp)
endforeach(prog)


add_subdirectory(validated_numerics)
add_subdirectory(optimisers)

