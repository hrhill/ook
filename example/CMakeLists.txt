set(PROGRAMS version)

foreach(ex ${PROGRAMS})
    add_executable(${ex} ${ex}.cpp)
endforeach(ex)


add_subdirectory(validated_numerics)
add_subdirectory(optimisers)

