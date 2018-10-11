cmake_minimum_required (VERSION 3.0)

#set your own eigen
#set(PATH_TO_EIGEN /Users/jianhangxu/Repos/eigen-git-mirror)

list (APPEND CMAKE_CXX_FLAGS "-std=c++11")		# c++11

project(mit)
# The version number.
set (MIT_VERSION_MAJOR 0)
set (MIT_VERSION_MINOR 0)

include_directories(${PATH_TO_EIGEN})
include_directories(${CMAKE_CURRENT_LIST_DIR}/include)
link_directories(${CMAKE_CURRENT_LIST_DIR}/lib)

configure_file (
  "${CMAKE_CURRENT_LIST_DIR}/MITConfig.h.in"
  "${CMAKE_CURRENT_LIST_DIR}/include/MITConfig.h"
  )

file(GLOB_RECURSE MIT_SRC_CPP ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp)

add_executable(${PROJECT_NAME}.x ${MIT_SRC_CPP} ${MIT_GLO_CPP})


#[[
cmake_minimum_required (VERSION 3.0)

list (APPEND CMAKE_CXX_FLAGS "-std=c++11")		# c++11


#add_library(Hello src/hello.cpp)

project(mit)
# The version number.
set (MIT_VERSION_MAJOR 0)
set (MIT_VERSION_MINOR 0)

include_directories(${CMAKE_CURRENT_LIST_DIR}/include)
#include_directories(${CMAKE_CURRENT_LIST_DIR}/include/global)

link_directories(${CMAKE_CURRENT_LIST_DIR}/lib)
#aux_source_directory(${CMAKE_CURRENT_LIST_DIR}/src ${hello_src})

configure_file (
  "${CMAKE_CURRENT_LIST_DIR}/MITConfig.h.in"
  "${CMAKE_CURRENT_LIST_DIR}/include/MITConfig.h"
  )

file(GLOB_RECURSE MIT_SRC_CPP ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp)
#file(GLOB MIT_GLO_CPP ${CMAKE_CURRENT_LIST_DIR}/src/global/*.cpp)

#add_library(global ${MIT_GLO_CPP})
#foreach(f ${MIT_SRC_CPP})
# message("${f}\n")
#endforeach(f)

#add_executable(${PROJECT_NAME}.x ${hello_src})
add_executable(${PROJECT_NAME}.x ${MIT_SRC_CPP} ${MIT_GLO_CPP})

#target_link_libraries(${PROJECT_NAME}.x global)
]]
