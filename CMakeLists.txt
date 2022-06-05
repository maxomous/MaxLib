cmake_minimum_required(VERSION 3.5)


set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(MaxLib
    DESCRIPTION "A Collection of Generally Handy Functions"
    VERSION 1.0
    LANGUAGES CXX
)

############################################################
# Create a library
############################################################

#Generate the shared library from the library sources
add_library(${PROJECT_NAME} STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/MaxLib.h
    ${CMAKE_CURRENT_SOURCE_DIR}/File/File.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Geom/Geom.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Devices/ADS1115/ads1115.cpp
)

############################################################
# Install
############################################################

# Library
install (TARGETS ${PROJECT_NAME}
    LIBRARY DESTINATION lib)

# Install Main Header File
INSTALL(FILES "${CMAKE_CURRENT_SOURCE_DIR}/MaxLib.h" DESTINATION include) # INSTALL(FILES ...) or install(DIRECTORY ...) 

# Build list of header files to install from other directorys. Root == "."
set(HEADER_DIRS 
    "Devices/ADS1115"
    "File" 
    "Geom" 
    "Vector" 
)
## Add Source Files from the other directories
foreach(DIR ${HEADER_DIRS})
    # Find all source files & append to list
    if(DIR STREQUAL ".")
        file(GLOB HEADER_FILES_IN_FOLDER *.h)
    else()
        file(GLOB HEADER_FILES_IN_FOLDER ${CMAKE_CURRENT_SOURCE_DIR}/${DIR}/*.h)
    endif()
    list(APPEND HEADER_FILES ${HEADER_FILES_IN_FOLDER})
endforeach()
# Install Header files
INSTALL(FILES ${HEADER_FILES} DESTINATION include/${PROJECT_NAME}) # INSTALL(FILES ...) or install(DIRECTORY ...) 


############################################################
# Uninstall
############################################################

if(NOT TARGET uninstall)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()

############################################################
# Options
############################################################

target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic) 
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)


