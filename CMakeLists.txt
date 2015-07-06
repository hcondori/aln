cmake_minimum_required (VERSION 3.0.0)
project (ALN)

#set(CMAKE_CXX_COMPILER clang++)

IF(WIN32)
  FIND_PACKAGE(OpenCL REQUIRED)
ENDIF(WIN32)

subdirs(src)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#add_definitions(-DOCL_SUPPORT)
