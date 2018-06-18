# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import division, unicode_literals

__author__ = "Travis W. Kemper, Scott Sides"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"


import os
from functools import wraps
import shutil

    
def tearDown_dir(func):
    @wraps(func)
    def wrapper_function(*args,**kwargs):
#        shutil.rmtree(os.path.join(TEST_DIR,'materials'))
#        shutil.rmtree(os.path.join(TEST_DIR,'storage'))
#        shutil.rmtree(os.path.join(TEST_DIR,'scripts'))
#        shutil.rmtree(os.path.join(TEST_DIR,'scratch'))
        shutil.rmtree('materials')
        shutil.rmtree('storage')
        shutil.rmtree('scripts')
        shutil.rmtree('scratch')
        return func(*args,**kwargs)
    return wrapper_function
