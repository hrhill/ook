set(TESTS test_line_search
          test_steepest_descent
          test_newton
          test_nonlinear_cg
          test_bfgs
          test_optimise
          test_mcstep
          test_functions
          test_factorisations
)

foreach(test ${TESTS})
    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                                  ${LAPACK_LIBRARIES}
                                  ${BLAS_LIBRARIES}
                                  ook)
    add_test(${test} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${test})
endforeach(test)

add_subdirectory(finite_differences)
