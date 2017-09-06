# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

import os
import sys
import platform

from setuptools import setup, find_packages, Extension

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='streamm',
    version='0.3.2',
    setup_requires=['numpy', 'setuptools>=18.0','pytest-runner'],    
    packages=find_packages(),
    # 
    url='http://streamm.nrel.gov',
    author='Travis W. Kemper, Ross Larsen and Scott W. Sides',
    author_email='travis.kemper.w@gmail.com',
    maintainer="Ross Larsen",
    maintainer_email="streamm@nrel.gov",
    #
    tests_require=['pytest'],
    
    # test_suite = ['streamm.structures.tests','streamm.forcefields.tests'],
    # data_files=['streamm/periodic_table.json']
    # packages=['streamm'],    
    # package_dir={'streamm': 'streamm'},
    # package_data={'streamm': ['periodic_table.json']},
    install_requires=['pymatgen'],
    ) 

