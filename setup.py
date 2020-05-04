#!/usr/bin/env python3
from setuptools import setup
from distutils.core import setup, Extension
import os
from external.get_include_paths import *
 
module_glassfire =\
    Extension(
        'glassfire', 
        sources = [
            'glassfire.cpp'
        ],
        include_dirs=[
            get_pybind_include(user=True),
            get_include(
                from_url='https://dl.bintray.com/boostorg/release/1.72.0/source/boost_1_72_0.tar.gz',
                to_local='external/_src/',
                extracted_folder='boost_1_72_0/',
                extract_method='gz'
            ),
            get_include(
                from_url='https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz',
                to_local='external/_src/',
                extracted_folder='eigen-3.3.7/',
                extract_method='gz'
            ),
        ],
        extra_compile_args=['-std=c++11', '-fpermissive'],
    )
 
setup (
    name = 'glassfire',
    version = '0.1.0',
    author = 'Yuan-Yen Tai',
    author_email = 'anrris330@gmail.com',
    ext_modules = [module_glassfire],
)
