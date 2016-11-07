file(GLOB_RECURSE SRC_FILES
     ${PROJECT_SOURCE_DIR}/test/*.cpp
     ${PROJECT_SOURCE_DIR}/src/*.cpp)

add_custom_target(cppcheck
                  COMMAND cppcheck
		  --enable=all
		  --std=c++11
		  --verbose
		  --quiet
		  -I${PROJECT_SOURCE_DIR}/include
		  ${SRC_FILES})
