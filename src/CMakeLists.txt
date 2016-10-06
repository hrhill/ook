add_library(ook SHARED
    vector.cpp
    matrix.cpp
)
target_link_libraries(ook ${BLAS_LIBRARIES})
