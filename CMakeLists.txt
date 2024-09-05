cmake_minimum_required(VERSION 3.6)
project(fastq_pair)

set(CMAKE_C_STANDARD 99)

# Add the zlib library from the external folder
add_subdirectory(external/zlib-1.3.1)

# List your source files
set(SOURCE_FILES main.c robstr.c fastq_pair.c is_gzipped.c is_gzipped.h)

# Add the executable for your project
add_executable(fastq_pair ${SOURCE_FILES})

# Link zlib (built locally) with your executable
target_link_libraries(fastq_pair PRIVATE zlibstatic)

# Installation configuration
install(TARGETS fastq_pair DESTINATION bin)
