#!/usr/bin/env python3
from setuptools import setup
from distutils.core import setup, Extension
import os
from web2local import *
from pathlib import Path
 
module_glassfire =\
    Extension(
        'glassfire', 
        sources = [ 'glassfire.cpp' ],
        include_dirs=[
            pybind11_path(user=True),
            #r'C:\\Users\\anrri\\Git\\Glassfire\\external\\boost_1_72_0\\boost_1_72_0\\',
            #r'C:\\Users\\anrri\\Git\\glassfire\\external\\eigen-3.3.7\\eigen-3.3.7\\',
            web2local(
                from_url='https://dl.bintray.com/boostorg/release/1.72.0/source/boost_1_72_0.tar.gz',
                to_local='external',
                extract_to='boost_1_72_0/',
                extract_method='gz'
            ),
            web2local(
                from_url='https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz',
                to_local='external',
                extract_to='eigen-3.3.7',
                extract_method='gz'
            )
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
