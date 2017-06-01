file(GLOB_RECURSE SRC_FILES ${PROJECT_SOURCE_DIR}/*.cpp)
file(GLOB_RECURSE HPP_FILES ${PROJECT_SOURCE_DIR}/*.hpp)

set(FORMAT_FILES ${HPP_FILES} ${H_FILES} ${SRC_FILES})

add_custom_target(
        clang_format
        COMMAND clang-format
        -i
        ${FORMAT_FILES}
)
