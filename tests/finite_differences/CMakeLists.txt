set(TESTS
   finite_differences 
)

include_directories(${CMAKE_SOURCE_DIR}/include/finite_differences)

foreach(test ${TESTS})
    add_executable(test_${test} test_${test}.cpp test_functions.cpp ${HEADERS})
    target_link_libraries(test_${test} ${Boost_TEST_EXEC_MONITOR_LIBRARY})
    add_test(test_${test} ${EXECUTABLE_OUTPUT_PATH}/test_${test})    
endforeach(test)

foreach(test ${TESTS})
    add_executable(parallel_test_${test} test_${test}.cpp test_functions.cpp ${HEADERS})
    target_link_libraries(parallel_test_${test} ${Boost_TEST_EXEC_MONITOR_LIBRARY})
    set_target_properties(parallel_test_${test} PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_PARALLEL ${OpenMP_CXX_FLAGS}" 
                                                            LINK_FLAGS "-D_GLIBCXX_PARALLEL ${OpenMP_CXX_FLAGS}")    
    add_test(parallel_test_${test} ${EXECUTABLE_OUTPUT_PATH}/parallel_test_${test})    
endforeach(test)

