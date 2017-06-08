# CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8)
project (fclib-cpp)

option(BUILD_DEMOS "Build test demos" ON)

SET(FCLIB_CPP_HEADERS
	fclib_class.hpp
)

SET(FCLIB_CPP_INCLUDEDIR
	${CMAKE_CURRENT_SOURCE_DIR}
)

find_package(HDF5 REQUIRED COMPONENTS HL)
include_directories(${HDF5_INCLUDE_DIRS})

set(FCLIB_CPP_LIBS ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

# Recurse into the "Hello" and "Demo" subdirectories. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
IF (BUILD_DEMOS)
	add_subdirectory (demo)
ENDIF()

