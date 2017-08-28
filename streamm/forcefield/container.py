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
parameters module

Class data structures force-field parameters

Most molecular dynamics codes set each particle to certain type (fftype)
usually each interaction is set using these fftypes

"""


class Container(object):
    """
    Container for force-field parameters
    """

    def __init__(self, verbose=False):
        """
        Constructor for Parameter.Container

        """
        self.ljtypes = dict()                           # Creates empty dict struc
        self.bondtypes = dict()                                   # Creates empty dict struc
        self.angletypes = dict()                                # Creates empty dict struc
        self.dihtypes = dict()                                # Creates empty dict struc
        self.imptypes = dict()                                  # Creates empty dict struc

        # Int count of the length of each dictionary
        #   mostly for internal use 
        self.n_ljtypes = 0    
        self.n_bondtypes = 0    
        self.n_angletypes = 0    
        self.n_dihtypes = 0    
        self.n_imptypes = 0            
        #
        # Set defaults 
        #
        self.nbfunc = 1      # Use 1 (Lennard-Jones) or 2 (Buckingham)
        self.combmixrule = 3   # 1 or 3  - geometric; 2 - arithmetic
        self.genpairs = "yes"  # generates 1-4 parameters that are not present in the pair list 
        self.fudgeLJ = 1.0   # the factor by which to multiply Lennard-Jones 1-4 interactions, default 1
        self.fudgeQQ = 1.0   # the factor by which to multiply electrostatic 1-4 interactions, default 1
        

    def __del__(self):
        """
        Destructor, clears object memory
        """
        # dictionary 
        del self.ljtypes
        del self.bondtypes
        del self.angletypes
        del self.dihtypes
        del self.imptypes
        # 
        del self.nbfunc
        del self.combmixrule
        del self.genpairs
        del self.fudgeLJ
        del self.fudgeQQ
        
    def __str__(self):
        """
        'Magic' method for printng contents
        """

        strucStr =  "\n"
        strucStr += sperator_line
        strucStr += "    Parameters \n"
        strucStr += sperator_line
        strucStr += "      LJ parameters %d \n"%(self.n_ljtypes)
        strucStr += "      Bond parameters %d \n"%(self.n_bondtypes)
        strucStr += "      Angle parameters %d \n"%(self.n_angletypes)
        strucStr += "      Dihedral parameters %d \n"%(self.n_dihtypes)
        strucStr += "      Imporper Dihedral parameters %d \n"%(self.n_imptypes)
        return strucStr

    def __iadd__(self, other ):
        """
        'Magic' method to implement the '+=' operator 
        
        """

        if isinstance(other, Container):
            for ljtkey_i, ljtype_i  in other.ljtypes.iteritems():
                self.add_LJtype(ljtype_i,deepcopy = True )
            for btkey_i,bondtype_i  in other.bondtypes.iteritems():
                self.add_bondtype(bondtype_i,deepcopy = True )
            for atkey_i,angletype_i  in other.angletypes.iteritems():
                self.add_angletype(angletype_i,deepcopy = True )
            for dtkey_i, dihtype_i  in other.dihtypes.iteritems():    
                self.add_dihtype(dihtype_i,deepcopy = True )
            for itkey_i, imptype_i  in other.imptypes.iteritems():    
                self.add_imptype(imptype_i,deepcopy = True )


        else:
            print "1st arg should be a Parameter.Container object"
            raise TypeError
        

        return self


        
    def add_LJtype(self, ljtype_i, deepcopy = True ):
        """
        Add 'Ljtype' object to ljtypes dict in this container and update n_ljtypes accordingly
        """
        if isinstance(ljtype_i, LJtype):
            self.n_ljtypes = len(self.ljtypes)
            if( deepcopy ):
                self.ljtypes[self.n_ljtypes] = copy.deepcopy(ljtype_i) # index 0 -> (N-1)
            else:
                self.ljtypes[self.n_ljtypes] = ljtype_i # index 0 -> (N-1)
                
            self.n_ljtypes = len(self.ljtypes)
        else:
            print "Attempting to add non-Ljtype type to container"
            raise TypeError




    def add_bondtype(self, bondtype_i, deepcopy = True ):
        """
        Add 'Bondtype' object to bondtypes dict in this container and update n_bondtypes accordingly
        """
        if isinstance(bondtype_i, Bondtype):
            self.n_bondtypes = len(self.bondtypes)
            if( deepcopy ):
                self.bondtypes[self.n_bondtypes] = copy.deepcopy(bondtype_i) # index 0 -> (N-1)
            else:
                self.bondtypes[self.n_bondtypes] = bondtype_i # index 0 -> (N-1)
                
            self.n_bondtypes = len(self.bondtypes)
        else:
            print "Attempting to add non-Bondtype type to container"
            raise TypeError


    def add_angletype(self, angletype_i, deepcopy = True ):
        """
        Add 'Angletype' object to angletypes dict in this container and update n_angletypes accordingly
        """
        if isinstance(angletype_i, Angletype):
            self.n_angletypes = len(self.angletypes)
            if( deepcopy ):
                self.angletypes[self.n_angletypes] = copy.deepcopy(angletype_i) # index 0 -> (N-1)
            else:
                self.angletypes[self.n_angletypes] = angletype_i # index 0 -> (N-1)
                
            self.n_angletypes = len(self.angletypes)
        else:
            print "Attempting to add non-Angletype type to container"
            raise TypeError

    def add_dihtype(self, dihtype_i, deepcopy = True ):
        """
        Add 'Dihtype' object to dihtypes dict in this container and update n_dihtypes accordingly
        """
        if isinstance(dihtype_i, Dihtype):
            self.n_dihtypes = len(self.dihtypes)
            if( deepcopy ):
                self.dihtypes[self.n_dihtypes] = copy.deepcopy(dihtype_i) # index 0 -> (N-1)
            else:
                self.dihtypes[self.n_dihtypes] = dihtype_i # index 0 -> (N-1)
                
            self.n_dihtypes = len(self.dihtypes)
        else:
            print "Attempting to add non-Dihtype type to container"
            raise TypeError


    def add_imptype(self, imptype_i, deepcopy = True ):
        """
        Add 'Imptype' object to imptypes dict in this container and update n_imptypes accordingly
        """
        if isinstance(imptype_i, Imptype):
            self.n_imptypes = len(self.imptypes)
            if( deepcopy ):
                self.imptypes[self.n_imptypes] = copy.deepcopy(imptype_i) # index 0 -> (N-1)
            else:
                self.imptypes[self.n_imptypes] = imptype_i # index 0 -> (N-1)
                
            self.n_imptypes = len(self.imptypes)
        else:
            print "Attempting to add non-Imptype type to container"
            raise TypeError
