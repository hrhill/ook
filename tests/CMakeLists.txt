set(TESTS optimise
          more_thuente)

foreach(test ${TESTS})
    add_executable(${test}_test ${test}_test.cpp)
    add_test(${test}_test ${EXECUTABLE_OUTPUT_PATH}/${test}_test)
    target_link_libraries(${test}_test ${Boost_TEST_EXEC_MONITOR_LIBRARY})
endforeach(test)
