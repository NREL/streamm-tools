#!/usr/bin/env python

from setuptools import setup

#  packages=['streamm', 'streamm.jobs', 'streamm.jobs', 'streamm.mpi', 'streamm.simulation', 'streamm.structure'],
# dep = ['os','sys','shutil','time','datetime','json','math','numpy','copy','scipy','random','json','optparse','logging','pickle']
#  dep openbabel ( py27-openbabel )

setup(name='streamm',
      version='0.3',
      url='http://streamm.nrel.gov',
      author='Travis W. Kemper and Scott W. Sides',
      author_email='travis.kemper@nrel.gov',
      test_suite = 'tests',
      # data_files=['streamm/periodic_table.json']
      packages=['streamm'],
      package_dir={'streamm': 'streamm'},
      package_data={'streamm': ['periodic_table.json']},
      #install_requires=['numpy','copy'],
      )

