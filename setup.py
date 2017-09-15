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
    version='0.3.1',
    setup_requires=['numpy', 'setuptools>=18.0','pytest-runner'],    
    packages=find_packages(),
    # 
    url='http://streamm.nrel.gov',
    author='Dr. Scott W. Sides, Dr. Travis W. Kemper, Dr. Ross E. Larsen and Dr. Peter Graf',
    author_email='streamm@nrel.gov',
    maintainer="Ross Larsen",
    maintainer_email="streamm@nrel.gov",
    license="Apache License, Version 2.0",
    description="The Simulation Toolkit for Renewable Energy and Advanced Materials Modeling (STREAMM)",
    long_description="",
    classifiers="",
    # 
    tests_require=['pytest'],
    # test_suite = ['streamm.structures.tests','streamm.forcefields.tests'],
    # data_files=['pymatgen_core/core/periodic_table.json'],
    # packages=['streamm'],    
    # package_dir={'streamm': 'streamm'},
    package_data={'pymatgen_core/core': ['pymatgen_core/core/periodic_table.json']},
    include_package_data=True,
    # install_requires=['pymatgen'],
    ) 

