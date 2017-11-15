# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import unicode_literals

__author__ = "Travis W. Kemper, Ph.D."
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"
__version__ = "0.3.4"
__email__ = "organicelectronics@nrel.gov"
__status__ = "Beta"




import shutil
import os
from functools import wraps
import glob

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_PATH =  os.path.join(TEST_DIR,'..','..','..','templates')
PRECISION = 4


def setUp_streamm(func):
    @wraps(func)
    def wrapper_function(*args,**kwargs):        
        os.chdir(TEST_DIR)
        return func(*args,**kwargs)
    return wrapper_function
          
def tearDown_streamm(func):
    @wraps(func)
    def wrapper_function(*args,**kwargs):
        os.chdir(TEST_DIR)
        for calc_dir_name in ['materials','storage','scripts','scratch']:
            dir_path = os.path.join(TEST_DIR,calc_dir_name)
            if( os.path.isdir(dir_path) ):
                shutil.rmtree(dir_path,ignore_errors=True)

        for sufix in ["*.json","*.log","*.xyz","*.csv","*.pkl"]:
            files = glob.glob(sufix)
            for fl_name in files:
                os.remove(os.path.join(TEST_DIR,fl_name))
                
        return func(*args,**kwargs)
    return wrapper_function
   
