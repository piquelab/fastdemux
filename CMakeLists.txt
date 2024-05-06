cmake_minimum_required(VERSION 3.10)
project(fastdemux)


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

# Set paths to htslib library and include directories
# htslib include and library paths (update these paths if necessary)
include_directories(/wsu/el7/groups/piquelab/samtools/1.19/include)
link_directories(/wsu/el7/groups/piquelab/samtools/1.19/lib)

# Include zlib headers
include_directories(${ZLIB_INCLUDE_DIRS})

# Define executables
add_executable(fastdemux fastdemux.c common.c)
add_executable(fastaseq fastaseq.c common.c)


# Link htslib
target_link_libraries(fastdemux PUBLIC hts m OpenMP::OpenMP_C ${ZLIB_LIBRARIES})
target_link_libraries(fastaseq PUBLIC hts m OpenMP::OpenMP_C ${ZLIB_LIBRARIES})


# Set Debug and Release configuration flags
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -pg -g")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")

set(CMAKE_INSTALL_PREFIX "/wsu/home/groups/piquelab/apps/el7/misc/")


# Set this before any add_executable or add_library calls
##set(CMAKE_SKIP_INSTALL_RPATH TRUE)
##set(CMAKE_INSTALL_RPATH "/wsu/el7/groups/piquelab/samtools/1.19/lib")

# Specify the install rules
install(TARGETS fastdemux fastaseq DESTINATION bin)
