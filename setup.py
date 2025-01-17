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
            web2local(
                from_url="https://boostorg.jfrog.io/artifactory/main/release/1.72.0/source/boost_1_72_0.tar.gz",
                to_local='.cpp_external',
                extract_to='boost_1_72_0',
                extract_method='gz',
                using_user_home=True
            ),
            web2local(
                from_url='https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz',
                to_local='.cpp_external',
                extract_to='eigen-3.3.7',
                extract_method='gz',
                using_user_home=True
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
