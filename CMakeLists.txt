#
# Created by Tai, Yuan-yen on 11/17/17.
# All rights reserved.
#

cmake_minimum_required(VERSION 3.8)
project(Glassfire)

set (CMAKE_CXX_STANDARD 11)

set(test_classifier test_classifier.cpp)

find_package(Boost 1.71.0)

set(INCLUDE_DIR
        /usr/local/Cellar/boost/1.71.0/include
        /usr/local/Cellar/eigen/3.3.7/include/eigen3
        )
include_directories(${INCLUDE_DIR})

add_executable(test_classifier ${test_classifier})

