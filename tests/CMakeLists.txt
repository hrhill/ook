set(TESTS optimise
          more_thuente
          test_functions
          factorisation)

foreach(test ${TESTS})
    add_executable(${test}_test ${test}_test.cpp)
    target_link_libraries(${test}_test ${Boost_TEST_EXEC_MONITOR_LIBRARY})
    add_test(${test} ${EXECUTABLE_OUTPUT_PATH}/${test}_test)    
endforeach(test)

add_subdirectory(finite_differences)
