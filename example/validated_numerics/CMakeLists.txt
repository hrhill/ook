add_executable(more_thuente94 more_thuente94.cpp)

add_executable(sdrive sdrive.cpp)
add_executable(sdrive_diagco sdrive_diagco.cpp)
target_link_libraries(sdrive ${OpenBLAS_LIBRARIES})
target_link_libraries(sdrive_diagco ${OpenBLAS_LIBRARIES})
