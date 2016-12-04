#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 16:38:53 2016

@author: sdebuyl
"""

from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
    ext_modules = cythonize('odetest.pyx'),  
    include_dirs=[numpy.get_include()]
)