#
# Created by Tai, Yuan-yen on 11/17/17.
# All rights reserved.
#

cmake_minimum_required(VERSION 3.8)
project(LocalGaussianModel)

set (CMAKE_CXX_STANDARD 11)


set(performance_test
        src/lgmbin/perf-test.cpp)

set(api_demo
        src/lgmbin/api-demo.cpp)

set(SOURCE_FILES
        src/cpp/lgm_classifier.h
        src/cpp/lgm_factory.h
        src/cpp/lgm_base.h
        include/lgm.h)

set (PYTHON_FILES
        src/py/cluster_generator.py
        src/py/plot_scatter.py)

set(BOOST_ROOT /usr/local/Cellar/boost/1.65.0)
find_package(Boost 1.65.0)

set(INCLUDE_DIR
        /usr/local/Cellar/eigen/3.3.4/include/eigen3
        include)
include_directories(${INCLUDE_DIR})

add_executable(perf-test ${performance_test} ${SOURCE_FILES} ${PYTHON_FILES})
add_executable(api-demo ${api_demo} ${SOURCE_FILES})

