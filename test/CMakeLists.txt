set(TESTS test_line_search
          test_steepest_descent
          test_fletcher_reeves
          test_bfgs
          test_optimise
          test_mcstep
          test_functions)

foreach(test ${TESTS})
    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} ${Boost_TEST_EXEC_MONITOR_LIBRARY} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    add_test(${test} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${test})
endforeach(test)

add_subdirectory(finite_differences)
