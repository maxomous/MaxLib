cmake_minimum_required(VERSION 3.5)


set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(File
    DESCRIPTION "Simplifies File IO Functions"
    VERSION 1.0
    LANGUAGES CXX
)

############################################################
# Create a library
############################################################

#Generate the shared library from the library sources
add_library(${PROJECT_NAME} STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/File.cpp
)

target_include_directories(${PROJECT_NAME}
    PUBLIC 
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic) 
