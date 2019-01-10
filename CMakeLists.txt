# Cmake test file 

cmake_minimum_required(VERSION 2.8.11)
project(RRR)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake-modules/")
#下面这个设置是自己加的，为了支持c++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fPIC -ggdb")


find_package(G2O)
find_package(Eigen3)
find_package(SuiteSparse)
find_package(Cholmod REQUIRED)
#find_package(gflags REQUIRED)
#add_subdirectory(gflags)

#set(GLOG_DIR /home/zihao/cprogram/code/packages/glog/rebuild/glog/)   #添加CaffeConfig.cmake的搜索路径

find_package(Glog REQUIRED)
find_package(GFlags REQUIRED)

if (NOT GLOG_FOUND)
    message(FATAL_ERROR "GLOG Not Found!")
endif (NOT GLOG_FOUND)
if ( GLOG_FOUND)
    include_directories(${GLOg_INCLUDE_DIR})
    message("GOOD news : GLOG is Found!")
    message("GLOG lib is found:" ${GLOG_INCLUDE_DIR})
endif ( GLOG_FOUND)

if (NOT GFlags_FOUND)
    message(FATAL_ERROR "GFlags Not Found!")
endif (NOT GFlags_FOUND)
if ( GFlags_FOUND)
    include_directories(${GFlags_INCLUDE_DIR})
    message("GOOD news : GFlags is Found!")
    message("GFlags lib is found:" ${GFlags_INCLUDE_DIR})
endif ( GFlags_FOUND)

#find_package(Glog REQUIRED)
#include_directories(BEFORE ${GLOG_INCLUDE_DIRS}) 

#find_package(glog)
#include_directories( /home/zihao/cprogram/code/packages/glog/rebuild/glog/ )

IF(CHOLMOD_FOUND)
    include_directories(${CHOLMOD_INCLUDE_DIR})
    message("CHOLMOD lib is found:" ${CHOLMOD_INCLUDE_DIR})
ENDIF(CHOLMOD_FOUND)

SET(CMAKE_BUILD_TYPE RelWithDebInfo)
#SET(CMAKE_BUILD_TYPE Debug)
if("${CMAKE_BUILD_TYPE}" STREQUAL "")
  #set(CMAKE_BUILD_TYPE Release)
  set(CMAKE_BUILD_TYPE Debug)
endif()


add_subdirectory(examples)
#add_subdirectory(glog)

