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

def readme():
    with open('description.txt') as f:
        return f.read()


setup(
    name='streamm',
    version='0.3.4',
    setup_requires=['setuptools>=18.0','pytest-runner',
                    "numpy>=1.9", "six", "requests", "ruamel.yaml>=0.15.6",
                      "monty>=0.9.6", "scipy>=0.14", "pydispatcher>=2.0.5",
                      "tabulate", "spglib>=1.9.9.44","pandoc","jupyter",
                      "matplotlib>=1.5", "palettable>=2.1.1", "sympy"],        
    packages=find_packages(),
    # 
    url='http://streamm.nrel.gov',
    author='Dr. Scott W. Sides, Dr. Travis W. Kemper, Dr. Ross E. Larsen and Dr. Peter Graf',
    author_email='organicelectronics@nrel.gov',
    maintainer="The National Renewable Energy Laboratory",
    maintainer_email="organicelectronics@nrel.gov",
    license="Apache License, Version 2.0",
    description="The Simulation Toolkit for Renewable Energy and Advanced Materials Modeling",
    long_description=readme(),
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Chemistry',
        ],
    # 
    tests_require=['pytest'],
    # test_suite = ['streamm.structures.tests','streamm.forcefields.tests'],
    # data_files=['pymatgen_core/core/periodic_table.json'],
    # packages=['streamm'],    
    # package_dir={'streamm': 'streamm'},
    package_data={'': ['templates/*.sh'], 
        '': ['templates/*.pbs'],
        '': ['templates/*.com'],
        '': ['templates/*.nw'],
        '': ['templates/*.in'],
        '': ['templates/*.template'],
        'pymatgen_core/core': ['pymatgen_core/core/periodic_table.json']},
    include_package_data=True,
    install_requires=["ruamel.yaml>=0.15.6","monty>=0.9.6","numpy>=1.9","scipy>=0.14"],
    keywords=[ "NWChem", "gaussian", "LAMMPS", "MD", "materials", "project","electronic", "structure", "analysis", "coupling", "molecular","dynamics", "organic", "generator"
                ,"high", "throughput", "combinatorial", "computational", "chemistry", "functionalization", "polymerization" ]
    ) 

