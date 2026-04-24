cmake_minimum_required(VERSION 3.10)
project(fastdemux)

include(GNUInstallDirs)


## For the debug build:
## cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug

# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()



# Set C standard
set(CMAKE_C_STANDARD 99)

# Find zlib and OpenMP
find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)

# Include directories for klib macro library
include_directories(${CMAKE_SOURCE_DIR}/klib)

# Discover htslib from the active environment (for example via loaded samtools module)
find_program(SAMTOOLS_EXECUTABLE samtools)
if(SAMTOOLS_EXECUTABLE)
  get_filename_component(SAMTOOLS_BIN_DIR "${SAMTOOLS_EXECUTABLE}" DIRECTORY)
  get_filename_component(SAMTOOLS_PREFIX "${SAMTOOLS_BIN_DIR}/.." ABSOLUTE)
endif()

find_path(HTSLIB_INCLUDE_DIR
  NAMES htslib/sam.h
  HINTS ${SAMTOOLS_PREFIX} $ENV{SAMTOOLS_ROOT} $ENV{SAMTOOLS_DIR}
  PATH_SUFFIXES include htslib/include
)

find_library(HTSLIB_LIBRARY
  NAMES hts libhts
  HINTS ${SAMTOOLS_PREFIX} $ENV{SAMTOOLS_ROOT} $ENV{SAMTOOLS_DIR}
  PATH_SUFFIXES lib lib64 htslib/lib
)

if(NOT HTSLIB_INCLUDE_DIR)
  message(FATAL_ERROR "Could not find htslib headers (missing htslib/sam.h). Set CMAKE_PREFIX_PATH or HTSLIB_INCLUDE_DIR.")
endif()

if(NOT HTSLIB_LIBRARY)
  message(FATAL_ERROR "Could not find libhts. Set CMAKE_PREFIX_PATH or HTSLIB_LIBRARY.")
endif()

get_filename_component(HTSLIB_LIBRARY_DIR "${HTSLIB_LIBRARY}" DIRECTORY)

# Include zlib headers
include_directories(${ZLIB_INCLUDE_DIRS})

# Define executables
add_executable(fastdemux fastdemux.c common.c)
add_executable(fastaseq fastaseq.c common.c)

target_include_directories(fastdemux PRIVATE ${HTSLIB_INCLUDE_DIR})
target_include_directories(fastaseq PRIVATE ${HTSLIB_INCLUDE_DIR})

set_target_properties(fastdemux PROPERTIES
  BUILD_RPATH "${HTSLIB_LIBRARY_DIR}"
  INSTALL_RPATH "${HTSLIB_LIBRARY_DIR}"
)

set_target_properties(fastaseq PROPERTIES
  BUILD_RPATH "${HTSLIB_LIBRARY_DIR}"
  INSTALL_RPATH "${HTSLIB_LIBRARY_DIR}"
)


# Link htslib
target_link_libraries(fastdemux PUBLIC ${HTSLIB_LIBRARY} m OpenMP::OpenMP_C ${ZLIB_LIBRARIES})
target_link_libraries(fastaseq PUBLIC ${HTSLIB_LIBRARY} m OpenMP::OpenMP_C ${ZLIB_LIBRARIES})


# Set Debug and Release configuration flags
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -pg -g")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")

# Keep installs local by default; callers can override with -DCMAKE_INSTALL_PREFIX
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install" CACHE PATH "Install path prefix" FORCE)
endif()


# Ensure installed binaries can find libhts without needing global env changes
set(CMAKE_INSTALL_RPATH "${HTSLIB_LIBRARY_DIR}")
set(CMAKE_BUILD_RPATH "${HTSLIB_LIBRARY_DIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Set this before any add_executable or add_library calls
##set(CMAKE_SKIP_INSTALL_RPATH TRUE)
##set(CMAKE_INSTALL_RPATH "/wsu/el7/groups/piquelab/samtools/1.19/lib")

# Specify the install rules
install(TARGETS fastdemux fastaseq DESTINATION ${CMAKE_INSTALL_BINDIR})
