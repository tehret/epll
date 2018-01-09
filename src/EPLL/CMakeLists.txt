
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
include_directories(${CMAKE_SOURCE_DIR}/src/Utilities)

add_library(epll SHARED
	epll.cpp
	LibMatrix.cpp
	../Utilities/io_png.c
	../Utilities/iio.c
	../Utilities/LibImages.cpp
	../Utilities/LibVideoT.cpp
	../Utilities/mt19937ar.c
	../Utilities/Utilities.cpp
	)


set_target_properties(epll PROPERTIES
	SOVERSION ${PROJECT_VERSION})

install(TARGETS epll
    LIBRARY DESTINATION lib
    COMPONENT libs)
