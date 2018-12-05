
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/EPLL)

add_executable(main_epll main_epll.cpp)
add_executable(comp_res compute_results.cpp)

message("-- LAPACK libs: " ${LAPACK_LIBRARIES})
message("-- LAPACK incs: " ${LAPACK_INCLUDES})
include_directories(${LAPACK_INCLUDES})
include_directories(${CBLAS_INCLUDE_DIR})

link_directories(${CMAKE_BINARY_DIR}/lib)
target_link_libraries(main_epll epll png ${LAPACK_LIBRARIES} ${CBLAS_LIBRARIES})
target_link_libraries(comp_res epll png ${LAPACK_LIBRARIES} ${CBLAS_LIBRARIES})
