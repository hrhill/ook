
cmake_minimum_required(VERSION 2.6)
project(OOK)

include(CheckIncludeFiles)
include(CheckLibraryExists)
include(CTest)
include(ExternalProject)

set(BUILD_SHARED_LIBS ON)

find_package(Boost 1.46 COMPONENTS test_exec_monitor program_options REQUIRED)
find_package(Doxygen)
find_package(OpenMP)

set(OpenBLAS_LIBRARIES libopenblas.so)
set(CMAKE_CXX_COMPILER "g++-5")

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Werror -Wall -pedantic -std=c++1y ${OpenMP_CXX_FLAGS}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
include_directories("${PROJECT_BINARY_DIR}/include")

add_subdirectory(example)
add_subdirectory(test)
add_subdirectory(src)

# ------------------------------------ Documentation setup --------------------
if(DOXYGEN_FOUND)
configure_file(${PROJECT_SOURCE_DIR}/Doxyfile.in
                   ${PROJECT_BINARY_DIR}/Doxyfile)
    add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM)

endif(DOXYGEN_FOUND)

#------------------------------------- Installation details -------------------
set(LIB_INST_DIR "${CMAKE_INSTALL_PREFIX}/lib64")
set(HEADER_INST_DIR "${CMAKE_INSTALL_PREFIX}/include")

install(DIRECTORY "${OOK_SOURCE_DIR}/include/"
    DESTINATION "${HEADER_INST_DIR}"
    FILES_MATCHING PATTERN "*.hpp"
    PATTERN "CMakeFiles" EXCLUDE
    PATTERN "Testing" EXCLUDE)

install(DIRECTORY "${PROJECT_BINARY_DIR}/lib"
        DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64}")
