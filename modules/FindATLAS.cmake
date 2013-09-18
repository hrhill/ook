FIND_PATH(ATLAS_INCLUDE_DIR cblas.h PATHS /usr/include /usr/local/include /usr/local/atlas/include)
FIND_PATH(ATLAS_LIB_DIR libtatlas.so PATHS /usr/lib64 
                                           /usr/local/lib64 
                                           /usr/lib 
                                           /usr/local/lib
                                           /usr/atlas/lib
                                           /usr/local/atlas/lib)

FIND_LIBRARY(ATLAS_LIBRARY NAMES tatlas PATHS ${ATLAS_LIB_DIR})

IF(ATLAS_INCLUDE_DIR AND ATLAS_LIBRARY)
   SET(ATLAS_FOUND TRUE)

ENDIF(ATLAS_INCLUDE_DIR AND ATLAS_LIBRARY)

IF(ATLAS_FOUND)
   IF (NOT ATLAS_FIND_QUIETLY)
      MESSAGE(STATUS "Found ATLAS: ${ATLAS_LIBRARY}")
   ENDIF (NOT ATLAS_FIND_QUIETLY)
   ADD_DEFINITIONS(-DHAVE_ATLAS)
   LINK_DIRECTORIES(${ATLAS_LIB_DIR})
   SET(ATLAS_LIBRARIES tatlas)
ELSE(ATLAS_FOUND)
   IF(ATLAS_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find ATLAS")
   ENDIF (ATLAS_FIND_REQUIRED)
   SET(ATLAS_LIBRARIES "")
ENDIF(ATLAS_FOUND)
