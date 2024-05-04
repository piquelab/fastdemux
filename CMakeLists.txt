cmake_minimum_required(VERSION 3.10)
project(fastdemux)

# Set C standard
set(CMAKE_C_STANDARD 99)

# Find zlib
find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)

# Include directories for htslib and klib
include_directories(${CMAKE_SOURCE_DIR}/klib)

# Set paths to htslib library and include directories
# htslib include and library paths (update these paths if necessary)
include_directories(/wsu/el7/groups/piquelab/samtools/1.19/include)
link_directories(/wsu/el7/groups/piquelab/samtools/1.19/lib)

# Include zlib headers
include_directories(${ZLIB_INCLUDE_DIRS})

# Your source files
add_executable(${PROJECT_NAME} main.c)

# Link htslib
target_link_libraries(${PROJECT_NAME} PUBLIC hts m OpenMP::OpenMP_C ${ZLIB_LIBRARIES})

# Set Debug and Release configuration flags
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -pg -g")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")

set(CMAKE_INSTALL_PREFIX "/wsu/home/groups/piquelab/apps/el7/misc/")

# Specify the install rules
install(TARGETS ${PROJECT_NAME} DESTINATION bin)
