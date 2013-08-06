
set(TESTS "")

foreach(test ${TESTS})
    add_executable(${test} EXCLUDE_FROM_ALL ${test}.cpp)
    target_link_libraries(test_${test}  ${Boost_TEST_EXEC_MONITOR_LIBRARY})
    add_test(${test} ${EXECUTABLE_OUTPUT_PATH}/${test})
endforeach(test)
