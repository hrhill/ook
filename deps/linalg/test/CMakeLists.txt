
set(TESTS norms
          linear_algebra
          factorisation
          operations)

set(LINK_LIBS ${LAPACK_LIBRARIES}
              ${BLAS_LIBRARIES}
              ${Boost_TEST_EXEC_MONITOR_LIBRARY})

foreach(test ${TESTS})
    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} ${LINK_LIBS})
    add_test(${test} ${EXECUTABLE_OUTPUT_PATH}/${test})
endforeach(test)
