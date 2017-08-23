# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

import os
import sys
import platform

from setuptools import setup, find_packages, Extension

#  packages=['streamm', 'streamm.jobs', 'streamm.jobs', 'streamm.mpi', 'streamm.simulation', 'streamm.structure'],
# dep = ['os','sys','shutil','time','datetime','json','math','numpy','copy','scipy','random','json','optparse','logging','pickle']
#  dep openbabel ( py27-openbabel )

setup(
    name='streamm',
    version='0.3.1',
    setup_requires=['numpy', 'setuptools>=18.0'],    
    packages=find_packages(),
    # 
    url='http://streamm.nrel.gov',
    author='Travis W. Kemper and Scott W. Sides',
    author_email='travis.kemper.w@gmail.com',
    maintainer="Ross Larsen",
    maintainer_email="streamm@nrel.gov",
    #
    test_suite = 'streamm.structure.tests',
    # data_files=['streamm/periodic_table.json']
    # packages=['streamm'],    
    # package_dir={'streamm': 'streamm'},
    # package_data={'streamm': ['periodic_table.json']},
    #install_requires=['numpy','copy'],
    ) 

