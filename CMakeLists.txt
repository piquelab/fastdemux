cmake_minimum_required(VERSION 3.10)
project(fastdemux)

# Set C standard
set(CMAKE_C_STANDARD 99)

find_package(OpenMP REQUIRED)
# Find zlib
find_package(ZLIB REQUIRED)


# Include directories for htslib and klib
include_directories(${CMAKE_SOURCE_DIR}/klib)

# Set paths to htslib library and include directories
# htslib include and library paths (update these paths if necessary)
include_directories(/wsu/el7/groups/piquelab/samtools/1.11/include)
link_directories(/wsu/el7/groups/piquelab/samtools/1.11/lib)

# Include zlib headers
include_directories(${ZLIB_INCLUDE_DIRS})

# Add htslib build
###add_subdirectory(htslib)

# Your source files
add_executable(${PROJECT_NAME} main.c)

# Link htslib
target_link_libraries(${PROJECT_NAME} PUBLIC hts OpenMP::OpenMP_C ${ZLIB_LIBRARIES})

