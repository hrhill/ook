set(TESTS
   test_finite_differences
)

foreach(test ${TESTS})
    add_executable(${test} ${test}.cpp test_functions.cpp)
    target_link_libraries(${test} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ook)
    add_test(${test} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${test})
endforeach(test)
