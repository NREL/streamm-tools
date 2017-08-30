.rsite# coding: utf-8
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
    # Import streamm dependstructures.encies
    from streamm.structures.containers import Container as strucCont
    # from streamm.structures.container import Container
    from streamm.structures.lattices import Lattice 
    from streamm.structures.nblists import NBlist 
    
    from streamm.structures.particle import Particle
    from streamm.structures.bonds import Bond 
    from streamm.structures.angles import Angle
    from streamm.structures.dihedrals import Dihedral
    from streamm.structures.impropers import Improper

except:
    print("streamm is not installed test will use relative path")
    import sys, os
    
    rel_path = os.path.join(os.path.dirname(__file__),'..')
    print("rel_path {}".format(rel_path))
    sys.path.append(rel_path)
    
    from bbatom import BBatom
 
    rel_path = os.path.join(os.path.dirname(__file__),'..','structure')
    print("rel_path {}".format(rel_path))
    sys.path.append(rel_path)
    
    from containers import Container as strucCont
    

class Attachment(object):
    '''
    Object to recode the attachment of two building block objects
    '''
    def __init__(self,tag_ij,tag_i,tag_j,rsite_i,n_i,rsite_j,n_j):

        self.tag_ij = tag_ij
        self.tag_i = tag_i
        self.tag_j = tag_j
        self.rsite_i = rsite_i
        self.n_i = n_i
        self.rsite_j = rsite_j
        self.n_j = n_j
    
    def __del__(self):
        """
        'Magic' method for deleting contents of container
        """
        del self.tag_ij
        del self.tag_i
        del self.tag_j
        del self.rsite_i
        del self.n_i
        del self.rsite_j
        del self.n_j

    def __str__(self):
        """
        'Magic' method for printng contents of container
        """
        return " %s ( %s %d ) +  %s ( %s %d ) -> %s "%(self.tag_i,self.rsite_i,self.n_i,self.tag_j,self.rsite_j,self.n_j,self.tag_ij)
    
class Container(strucCont):
    """
    Data structure for describing a section of a molecule
       that can be joined together to make a new molecule 
    """    
    def __init__(self,tag="blank"):    
        """
        Constructor for a composite structures. 
        """
        strucCont.__init__(self)
        #super(structures.Container,self).__init__()
        self.tag = tag 
        self.common_tag = tag 
        # Reactive sites 
        self.n_term = int(0)
        self.n_func = int(0)
        self.n_sub = int(0)
        self.terms = []
        self.funcs = []
        self.subs = []        
        self.attachments = []
        
    def __del__(self):
        del self.n_term
        del self.n_func
        del self.n_sub
        del self.terms
        del self.funcs
        del self.subs
        del self.attachments
        
    def calc_attachments(self):
        '''
        Count number terminal and functional attachment points
        '''
        #
        self.n_term = 0 
        self.n_func = 0
        self.n_sub = 0
        #
        for pkey_i, particle_i  in self.particles.iteritems():
            if( particle_i.rsite == "T" ):
                self.n_term += 1 
            elif( particle_i.rsite == "F" ):
                self.n_func += 1 
            elif( particle_i.rsite == "S" ):
                self.n_sub += 1
        #       
        logger.info("Found {} terminal sites, {} functionalizable sites and {} substitutional sites ".format(self.n_term,self.n_func,self.n_sub ))
        #
        return 
 

    def find_rsite(self,rsite_i,n_i,Xn_i):
        '''
        Find keys of terminal atom X and it's cap R
        
        Arguments:
            rsite_i   (str) Reactive site name 
            n_i       (int) number of connection to used in  container 1 
            Xn_i      (int) nieghbor number of R atom

        Return:
            Rkey  - the particle to be replaced
            XKey  - the nieghbor Xn_i of the particle to be replaced 
        '''
        cnt_i = 0

        # print " find_XR self.tag ",self.tag
        
        for pkey_i, particle_i  in self.particles.iteritems():
            
            if( particle_i.rsite == rsite_i ):
                # print " debug self.bonded_nblist.getnbs(pkey_i) ",cnt_i,self.bonded_nblist.getnbs(pkey_i)[Xn_i]
                if( cnt_i == n_i ):
                    Rkey_i = pkey_i
                    Xkey_i = self.bonded_nblist.getnbs(pkey_i)[Xn_i]
                    # R - X - StrucC 
                    return Rkey_i,Xkey_i
                cnt_i += 1
        error_line = "Particle key n=%d with bbid = %s and neighbor %d in building block %s not found "%(n_i,rsite_i,Xn_i,self.tag)
        raise LookupError(error_line)

    def align_yaxis(self, key_k,yangle):
        '''
        Align position of particle k along the y-axis

                         y
                         ^
                         |      key_k (First non H nieghbor of i)
                         |    / 
        x <---  R    - key_i

        yangle (float) angle in radians
   
        '''
        # Get position of i and j 
        pos_k = self.positions[key_k]

        logger.debug(".align_yaxis pos_k {} ".format(str(pos_k)))
        
        # Get the zy component of position k
        pos_yz = np.array([0.0,pos_k[1],pos_k[2]])
        # Unit vector in the y direction 
        unit_y = np.array( [0.0,1.0,0.0])

        # Rotate structure around x
        # according to angle between position k and y-axis
        # Subtract desired final  angle with the y-axis
        logger.debug("yangle:{}".format(yangle))
        cos_yz = calc_cos(unit_y,pos_yz)
        angle_yz = np.arccos(cos_yz)
        
        logger.debug(" angle_yz:{}".format(angle_yz))
        
        angle_yz += -1.0*yangle
        
        logger.debug("  angle_yz - yangle: {}".format(angle_yz))
        
        if( pos_yz[2] > 0.0 ):
            # if in the 1st or 2nd quadrents rotate  clockwise
            direction="clockwise"
        elif( pos_yz[2] <= 0.0 ):
            # if in the 3rd or 4th quadrents rotate counter clockwise 
            # of if 180 rotation
            direction="counterclockwise"
        logger.debug("    direction:{}".format(direction))
        self.rotate_yz(angle_yz,direction=direction,verbose=False)
        
        return

    def align_bond(self, key_i, key_j):
        '''
        Align bond along the x-axis

                         y
                         ^
                         |    
                         |     
        x <--- key_j - key_i - Struc
                              
   
        '''
        # Get position of i and j 
        pos_i = self.positions[key_i]
        pos_j = self.positions[key_j]

        logger.debug(".align_bond pos_i {} pos_j {}".format(pos_i,pos_j))

        # Shift particle i to origin 
        self.shift_pos(-1.0*pos_i)
        # Get the xy component of position j
        pos_xy = np.array([pos_j[0],pos_j[1],0.0])
        # Unit vector in the x direction 
        unit_x = np.array( [1.0,0.0,0.0])

        # Rotate structure around z
        # according to angle between position j and x-axis
        if( pos_xy[1] > 0.0 ):
            # if in the 1st or 2nd quadrents rotate  clockwise 
            self.rotate_xy(np.arccos(calc_cos(unit_x,pos_xy)),direction="clockwise",verbose=False)
        elif( pos_xy[1] < 0.0 ):
            # if in the 3rd or 4th quadrents rotate counter clockwise 
            # of if 180 rotation 
            self.rotate_xy(np.arccos(calc_cos(unit_x,pos_xy)),direction="counterclockwise",verbose=False)
        else:
            logger.debug(" Bond already aligned skipping rotation")
        # if  pos_xy[1] == 0.0 no rotation is needed

        #if( debug ):
        #print "     >align_bond after xy rot self.positions[key_j] ",self.positions[key_j]
        # Update position j after rotation 
        pos_j = np.array(self.positions[key_j]  )
        # Get components of position j 
        pos_xz = np.array([pos_j[0],0.0,pos_j[2]])
        #
        logger.debug(" pos_xz:{}".format(pos_xz))
        
        # Rotate structure around y
        # according to angle between position j and x-axis
        if( pos_xz[2] > 0.0 ):
            # if in the 1st or 2nd quadrents rotate  clockwise 
            self.rotate_xz(np.arccos(calc_cos(unit_x,pos_xz)),direction="clockwise",verbose=False)
        elif( pos_xz[2] <= 0.0 ):
            # if in the 3rd or 4th quadrents rotate counter clockwise 
            # of if 180 rotation 
            self.rotate_xz(np.arccos(calc_cos(unit_x,pos_xz)),direction="counterclockwise",verbose=False)
        # if  pos_xy[1] == 0.0 no rotation is needed

        logger.debug("After xz rot self.positions[key_j]:{}".format(self.positions[key_j] ))

        return
    
    def getSubStructure(self,pkeys,tag="blank"):
        """
        This is redundant to getSubStructure in structures.py
        """
        new_strucC = Container(str(tag))
        
        key_update = dict()
        # Set lattice 
        new_strucC.lat = self.lat
        # Set particles
        for pkey_i in pkeys:
            p_i = self.particles[pkey_i]
            pos_i = self.positions[pkey_i]            
            new_strucC.add_partpos(p_i,pos_i, deepcopy = True)
            key_update[pkey_i]  = new_strucC.n_particles -1
            
        if( len(self.bonded_nblist.index) > 0 ):
            # Update bonded nieghbor list
            new_strucC.bonded_nblist = structures.NBlist() 
            for pkey_i in pkeys:
                new_strucC.bonded_nblist.index.append(new_strucC.bonded_nblist.cnt + 1)
                for pkey_j in self.bonded_nblist.getnbs(pkey_i):
                    if( pkey_j in pkeys ):
                        new_strucC.bonded_nblist.cnt += 1 
                        new_strucC.bonded_nblist.list.append(key_update[pkey_j])

            new_strucC.bonded_nblist.index.append(new_strucC.bonded_nblist.cnt + 1)
            new_strucC.bonded_bonds()

        return new_strucC
    
    def prepattach(self,rsite_i="R",n_i=0,Xn_i=0,dir=1,yangle=0.0):
        '''
        Preppair a building block to be attached to another

        This is handy for batch attachments to not repeat the intial steps in the attach function

        Arguments:
            yangle (float) angle in radians
            
        '''
    
        logger.debug(".prepattach {} {} {}".format(rsite_i,n_i,Xn_i))
        
        if( dir != 1 and dir != -1 ):
            error_line = " dir is not 1 or -1 "
            raise  ValueError(error_line)
            
        
        bb_prepped = copy.deepcopy(self)
        if( bb_prepped.n_particles > 1 ):
            
            # Find keys of attachment points
            #   j - i - Struc
            # Rkey_i,Xkey_i
            key_j,key_i = bb_prepped.find_XR(rsite_i,n_i,Xn_i)
            # Align building blocks along bonds of attachment atoms
            if( dir == 1 ):
                bb_prepped.align_bond(key_i,key_j,debug = debug)
                pos_i = bb_prepped.positions[key_i]
                logger.debug("Direction:{} align_bond key:{} pos_i:{} ".format(dir,key_i,str(bb_prepped.positions[key_i])))
            elif( dir == -1 ):
                bb_prepped.align_bond(key_j,key_i,debug = debug)
                # 
                pos_i = bb_prepped.positions[key_i]
                bb_prepped.shift_pos(-1.0*pos_i)
                logger.debug("Direction:{} align_bond key:{} pos_i:{} ".format(dir,key_i,str(bb_prepped.positions[key_i])))

            # Find attached heavy atoms
            # align fist heavy neighbor with y axis
            for key_k in bb_prepped.bonded_nblist.getnbs(key_i):
                    particle_k = bb_prepped.particles[key_k]
                    if( particle_k.number != 1 ):                        
                            bb_prepped.align_yaxis(key_k,yangle,debug = debug)
                            break

            # Remove atoms at R in building block i
            bb_prepped.del_particle(key_j)
            Xo_i = bb_prepped.keyupdate[key_i]
            # record attachment point 
            bb_prepped.attach_p = Xo_i
        else:
            # This is a particle substitution
            pos_i = bb_prepped.positions[0]
            bb_prepped.shift_pos(-1.0*pos_i)
            bb_prepped.attach_p = 0
            

        return bb_prepped

