set(TESTS line_search
          optimise
          test_functions)

foreach(test ${TESTS})
    add_executable(${test}_test ${test}.cpp)
    target_link_libraries(${test}_test ${Boost_TEST_EXEC_MONITOR_LIBRARY} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    add_test(${test} ${EXECUTABLE_OUTPUT_PATH}/${test}_test)    
endforeach(test)

add_subdirectory(finite_differences)
