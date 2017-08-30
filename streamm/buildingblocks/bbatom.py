# coding: utf-8
# Copyright (c) Alliance for Sustainable Energy, LLC
# Distributed under the terms of the Apache License, Version 2.0

from __future__ import division, unicode_literals

__author__ = "Travis W. Kemper, Scott Sides"
__copyright__ = "Copyright 2015, Alliance for Sustainable Energy, LLC"
__version__ = "0.3"
__email__ = "streamm@nrel.gov"
__status__ = "Beta"

"""
This module defines the classes relating to a molecular building block 
"""

import numpy as np 
import copy
import sys



import logging
logger = logging.getLogger(__name__)


try:
    
    from streamm.structure.atoms import Atom

except:
    print("streamm is not installed test will use relative path")
    import sys, os
    rel_path = os.path.join(os.path.dirname(__file__),'..','structure')
    print("rel_path {}".format(rel_path))
    sys.path.append(rel_path)
    from atoms import Atom

logger.warning("The .cplytag attribute will be deprecated in v0.3.2 ")
        
class BBatom(Atom):
    """
    A derived type of particles for atoms in a building block

    bb - X - T  - terminal site for polymerization 
    bb - X - R  - functionalizable site for functional groups 
    bb - X - S  - substitutable site for atomic substitutions 

    """

    def __init__(self,type="bbatom",symbol="X"):
        Atom.__init__(self,type=type, symbol=symbol)
        self.bbid= ""
        self.cplytag =  ""
        

    def __del__(self):
        """
        Destructor, clears structure memory and calls held container destructors
        """
        del self.bbid
        del self.cplytag
            

    def __str__(self):
        """
        'Magic' method for printng contents of container
        """
        return " %s (%s)"%(self.tag,self.bbid )
