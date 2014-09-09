
cmake_minimum_required(VERSION 2.6)
project(OOK)

include(CheckIncludeFiles)
include(CheckLibraryExists)
include(CTest)
include(ExternalProject)

set(BUILD_SHARED_LIBS ON)
set(CMAKE_CXX_FLAGS  "-Wall -pedantic -std=c++1y -D__STRICT_ANSI__")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

find_package(Boost 1.46 COMPONENTS test_exec_monitor program_options REQUIRED)
find_package(Doxygen)
find_package(OpenMP)
find_package(BLAS)
find_package(LAPACK)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
include_directories("${PROJECT_BINARY_DIR}/include")

add_subdirectory(example)
add_subdirectory(test)

#ExternalProject_Add(numeric_bindings
#    GIT_REPOSITORY "https://github.com/hrhill/numeric_bindings.git"
#    PREFIX "deps/numeric_bindings"
#    BINARY_DIR ${PROJ_BINARY_DIR}
#    CMAKE_ARGS
#      -DCMAKE_INSTALL_PREFIX:PATH=${PROJECT_BINARY_DIR}
#)

#ExternalProject_Add(linalg
#    GIT_REPOSITORY "https://github.com/hrhill/linalg.git"
#    GIT_TAG "v0.1.3"
#    PREFIX "deps/linalg"
#    BINARY_DIR ${PROJ_BINARY_DIR}
#    CMAKE_ARGS
#      -DCMAKE_INSTALL_PREFIX:PATH=${PROJECT_BINARY_DIR}
#      -DBLAS_LIBRARIES=${BLAS_LIBRARIES}
#      -DLAPACK_LIBRARIES=${LAPACK_LIBRARIES}
#    DEPENDS numeric_bindings
#)


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
