
cmake_minimum_required(VERSION 2.6)
project(OOK)

include(CheckIncludeFiles)
include(CheckLibraryExists)
include(CTest)

set(BUILD_SHARED_LIBS ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${OOK_SOURCE_DIR}/modules/")
set(CMAKE_CXX_FLAGS  "-Wall -pedantic -std=gnu++0x -D__STRICT_ANSI__")
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

find_package(Boost 1.50 COMPONENTS test_exec_monitor program_options REQUIRED)
find_package(Doxygen)
find_package(OpenMP)
find_package(BLAS)
find_package(LAPACK)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
add_subdirectory(example)
add_subdirectory(test)

# ------------------------------------ Documentation setup ------------------------------------------
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile)
    add_custom_target(doc ${DOXYGEN_EXECUTABLE}
                          ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
                          WORKING_DIRECTORY
                          ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif(DOXYGEN_FOUND)

#------------------------------------- Installation details -----------------------------------------
set(LIB_INST_DIR "${CMAKE_INSTALL_PREFIX}/lib64")
set(HEADER_INST_DIR "${CMAKE_INSTALL_PREFIX}/include")

install(DIRECTORY "${OOK_SOURCE_DIR}/include/"
    DESTINATION "${HEADER_INST_DIR}"
    FILES_MATCHING PATTERN "*.h"
    PATTERN "CMakeFiles" EXCLUDE
    PATTERN "Testing" EXCLUDE)
