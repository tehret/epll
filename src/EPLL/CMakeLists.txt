
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

add_library(epll SHARED
	epll.cpp
	LibMatrix.cpp
	iio.c
	LibImages.cpp
	mt19937ar.c
	)


set_target_properties(epll PROPERTIES
	SOVERSION ${PROJECT_VERSION})

install(TARGETS epll
    LIBRARY DESTINATION lib
    COMPONENT libs)
