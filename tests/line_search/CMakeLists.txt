set(TESTS more_thuente
          )

foreach(test ${TESTS})
    add_executable(${test}_test ${test}.cpp)
    target_link_libraries(${test}_test ${Boost_TEST_EXEC_MONITOR_LIBRARY})
    add_test(${test} ${EXECUTABLE_OUTPUT_PATH}/${test}_test)    
endforeach(test)

