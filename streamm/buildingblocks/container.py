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
This module defines the classes relating to a molecular building blocks 
"""

import numpy as np 
import copy
import sys



import logging
logger = logging.getLogger(__name__)


try:
    # Import streamm dependstructures.encies
    from streamm.structures.container import Container as strucCont
    # from streamm.structures.container import Container
    from streamm.structures.lattice import Lattice 
    from streamm.structures.nblist import NBlist 
    
    from streamm.structures.particle import Particle
    from streamm.structures.bond import Bond 
    from streamm.structures.angle import Angle
    from streamm.structures.dihedral import Dihedral
    from streamm.structures.improper import Improper

except:
    print("streamm is not installed test will use relative path")
    import sys, os
    
    rel_path = os.path.join(os.path.dirname(__file__),'..')
    print("rel_path {}".format(rel_path))
    sys.path.append(rel_path)
    
    # from bbatom import BBatom
 
    rel_path = os.path.join(os.path.dirname(__file__),'..','structure')
    print("rel_path {}".format(rel_path))
    sys.path.append(rel_path)
    
    from container import Container as strucCont
    from bond import Bond 
    from nblist import NBlist 


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
    def __init__(self,tag=str("blank"),matrix=[100.0,0.0,0.0,0.0,100.0,0.0,0.0,0.0,100.0]):
        """
        Constructor for a composite structures. 
        """
        strucCont.__init__(self, tag=tag, matrix=matrix)
        # Reactive sites 
        self.n_func = int(0)
        self.funcs = {}   
        self.attachments = []
        
    def __del__(self):
        del self.n_func
        del self.funcs
        del self.attachments
        # 
    def find_rsites(self):
        '''
        Find dictionary of lists of particle indexes based on rsites type
        '''
        #
        self.funcs = {}
        #
        for pkey_i, particle_i  in self.particles.iteritems():
            if( len(particle_i.rsite) > 0 ):
                if( particle_i.rsite not in self.funcs.keys() ):
                    self.funcs[particle_i.rsite] = []
                self.funcs[particle_i.rsite].append(pkey_i)
        #
        self.n_func = 0
        for type,list in self.funcs.iteritems():
            logger.info("Found  {} functionalizable sites with type {} ".format(len(list),type))
            self.n_func += len(list)
        #       
        logger.info("Found  {} functionalizable sites ".format(self.n_func))
        #
        return
    
    def show_rsites(self):
        '''
        Show the functional sites 
        '''
        out_line =''
        for rsite_i,rsite_list in self.funcs.iteritems():
            for pkey_i in rsite_list:
                out_line += "rsite:{}[ paticle:{} index:{} n_bonds:{}] \n".format(rsite_i,str(self.particles[pkey_i]),pkey_i,self.bonded_nblist.calc_nnab(pkey_i))
                
        return out_line 

    def get_rsite(self,rsite_i,n_i=0,Xn_i=0):
        '''
        Get particle key of reactive site and the paricle key of it's Xn_i connected neighbor 
        
        Arguments:
            rsite_i   (str) Reactive site type 
            n_i       (int) index of reactive site within in list of 
            Xn_i      (int) nieghbor number of reactive site 

        Return:
            Rkey  - the particle to be replaced
            XKey  - the nieghbor Xn_i of the particle to be replaced
            
        '''
        cnt_i = 0

        # print " get_rsite self.tag ",self.tag
        if( rsite_i in self.funcs.keys() ):
            
            rsite_list =  self.funcs[rsite_i]
            if( n_i < len(rsite_list) ):
                Rkey_i = rsite_list[n_i]
                Xkey_i = self.bonded_nblist.getnbs(Rkey_i)[Xn_i]
                return Rkey_i,Xkey_i
            else:
                raise KeyError("rsite index {} not found in list funcs[{}]  ".format(n_i,rsite_i))
                return None,None            
        else:
            raise KeyError("rsite {} not found in funcs dictionary make sure to run .find_rsites to polulate funcs ".format(rsite_i))
            return None,None
        

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
        self.rotate_yz(angle_yz,direction=direction)
        
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
            self.rotate_xy(np.arccos(calc_cos(unit_x,pos_xy)),direction="clockwise")
        elif( pos_xy[1] < 0.0 ):
            # if in the 3rd or 4th quadrents rotate counter clockwise 
            # of if 180 rotation 
            self.rotate_xy(np.arccos(calc_cos(unit_x,pos_xy)),direction="counterclockwise")
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
            self.rotate_xz(np.arccos(calc_cos(unit_x,pos_xz)),direction="clockwise")
        elif( pos_xz[2] <= 0.0 ):
            # if in the 3rd or 4th quadrents rotate counter clockwise 
            # of if 180 rotation 
            self.rotate_xz(np.arccos(calc_cos(unit_x,pos_xz)),direction="counterclockwise")
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
    
    def prepattach(self,rsite_i,n_i=0,Xn_i=0,dir=1,yangle=0.0):
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
            key_j,key_i = bb_prepped.get_rsite(rsite_i,n_i,Xn_i)
            # Align building blocks along bonds of attachment atoms
            if( dir == 1 ):
                bb_prepped.align_bond(key_i,key_j)
                pos_i = bb_prepped.positions[key_i]
                logger.debug("Direction:{} align_bond key:{} pos_i:{} ".format(dir,key_i,str(bb_prepped.positions[key_i])))
            elif( dir == -1 ):
                bb_prepped.align_bond(key_j,key_i)
                # 
                pos_i = bb_prepped.positions[key_i]
                bb_prepped.shift_pos(-1.0*pos_i)
                logger.debug("Direction:{} align_bond key:{} pos_i:{} ".format(dir,key_i,str(bb_prepped.positions[key_i])))

            # Find attached heavy atoms
            # align fist heavy neighbor with y axis
            for key_k in bb_prepped.bonded_nblist.getnbs(key_i):
                    particle_k = bb_prepped.particles[key_k]
                    if( particle_k.symbol != 'H' ):                        
                            bb_prepped.align_yaxis(key_k,yangle)
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


def calc_cos(v_i,v_j):
    """
    Calculate the cos theta between two vectors 
    """
    return np.dot(v_i/np.linalg.norm(v_i),v_j/np.linalg.norm(v_j))

def checkprep(bbC_i,bbC_j,covbuffer=1.5):
    '''
    Check 
    
    Arguments:
            bblockC_i (Container) Buildingblock container 1
            Xo_i (int) key of attachment point particle container 1
            bblockC_j (Container) Buildingblock container 2 
            Xo_j (int) key of attachment point particle container 2
    Retrun
        True - if no overlap
        False - if overlap is found
        
        
    '''
    #
    # NoteTK This seem backwards
    #
    # Set attachment points 
    Xo_i = bbC_i.attach_p
    Xo_j = bbC_j.attach_p
    
    npos_i = bbC_i.positions
    #npos_i.pop(Xo_i)
    npos_j = bbC_j.positions 
    #npos_j.pop(Xo_j)
    lat_i = bbC_i.lat
    npos_ij,nd_ij = lat_i.delta_npos(npos_i,npos_j)

    for pkey_i, particle_i  in bbC_i.particles.iteritems():
        if( pkey_i != Xo_i and particle_i.symbol != 'H' ):
            for pkey_j, particle_j  in bbC_j.particles.iteritems():
                if( pkey_j != Xo_j and particle_j.symbol != 'H' ):
                    dij = nd_ij[pkey_i][pkey_j]
                    radii_i = bbC_i.particles[pkey_i].bonded_radius
                    radii_j = bbC_j.particles[pkey_j].bonded_radius
                    cut_off = radii_i +radii_j
                    cut_off = cut_off*covbuffer
                    logger.debug(" i:{}-j:{} dr_ij:{} with cut off:{} ".format(pkey_i,pkey_j,dij,cut_off))
                    if( dij < cut_off ):
                        log_line = "      >checkprep \n"
                        log_line += "           particle i %s %d \n"%(particle_i.properties['symbol'],pkey_i)
                        log_line += "           particle j %s %d \n"%(particle_j.properties['symbol'],pkey_j)
                        log_line += "           radii_i %f \n"%(radii_i)
                        log_line += "           radii_j %f \n"%(radii_j)
                        log_line += "           dij %f \n"%(dij)
                        logger.debug(log_line)
                        return False 
                    
    return True 
    
def shiftprep(bblockC_i,bblockC_j,debug=False):
    '''

    Concatenate two building block containers that had prepattach already ran on them 

        This is handy for batch attachments to not repeat the intial steps in the attach function
    
    Arguments:
            bblockC_i (Container) Buildingblock container 1
            Xo_i (int) key of attachment point particle container 1
            bblockC_j (Container) Buildingblock container 2 
            Xo_j (int) key of attachment point particle container 2
            bbid_i   (str) Connecting atom bbid in container 1 
            n_i      (int) number of connection to used in  container 1 
            bbid_j   (str) Connecting atom bbid in container 2
            n_j      (int) number of connection to used in  container 2

        { bblockC_i - X_i - R_i } +  { R_j - X_j - bblockC_j }
                          =>
        { bblockC_i - X_i - X_j - bblockC_j }
    '''

    bbC_i = copy.deepcopy(bblockC_i)
    bbC_j = copy.deepcopy(bblockC_j)
    #
    Xo_i = bbC_i.attach_p
    Xo_j = bbC_j.attach_p
    # Shift  building block j to correct bond length 
    radii_i = bbC_i.particles[Xo_i].bonded_radius
    radii_j = bbC_j.particles[Xo_j].bonded_radius
    bond_vec = np.array([radii_i + radii_j,0.0,0.0])

    if( debug ):
        print " >attachprep Xo_i ",Xo_i
        print " >attachprep bbC_i.position[Xo_i] ", bbC_i.positions[Xo_i]
        print " >attachprep radii_i ",radii_i
        print " >attachprep Xo_j ",Xo_j
        print " >attachprep bbC_j.position[Xo_j] ", bbC_j.positions[Xo_j]
        print " >attachprep radii_j ",radii_j
    
    bbC_j.shift_pos(-1.0*bond_vec)

    return bbC_i,bbC_j
    
def attachprep(bbC_i,bbC_j):
    '''

    Concatenate two building block containers that had prepattach already ran on them 

        This is handy for batch attachments to not repeat the intial steps in the attach function
    
    Arguments:
            bblockC_i (Container) Buildingblock container 1
            Xo_i (int) key of attachment point particle container 1
            bblockC_j (Container) Buildingblock container 2 
            Xo_j (int) key of attachment point particle container 2
            bbid_i   (str) Connecting atom bbid in container 1 
            n_i      (int) number of connection to used in  container 1 
            bbid_j   (str) Connecting atom bbid in container 2
            n_j      (int) number of connection to used in  container 2

        { bblockC_i - X_i - R_i } +  { R_j - X_j - bblockC_j }
                          =>
        { bblockC_i - X_i - X_j - bblockC_j }
    '''
    # Set attachment points 
    Xo_i = bbC_i.attach_p
    Xo_j = bbC_j.attach_p
    #
    # Add j to i
    bbC_i += bbC_j
    #
    # Get updated atom key to form bond 
    Xp_j = bbC_i.keyupdate[Xo_j]
    #
    # Create bond between  X_i - X_j
    bond_ij = Bond(Xo_i,Xp_j)
    bbC_i.add_bond(bond_ij)
    #
    # Remake neighbor list based on updated bonds 
    bbC_i.bonded_nblist = NBlist()
    bbC_i.bonded_nblist = bbC_i.guess_nblist(0,radii_buffer=1.25)
    # Update number of attachment points
    bbC_i.find_rsites()                

    return bbC_i
    
        
def attach(bblockC_i,bblockC_j,bbid_i="R",n_i=0,bbid_j="R",n_j=0,tag="blank"):
        '''

        Concatenate two building block containers

        Arguments:
            bblockC_i (Container) Buildingblock container 1 
            bblockC_j (Container) Buildingblock container 2 
            bbid_i   (str) Connecting atom bbid in container 1 
            n_i      (int) number of connection to used in  container 1 
            bbid_j   (str) Connecting atom bbid in container 2
            n_j      (int) number of connection to used in  container 2

        { bblockC_i - X_i - R_i } +  { R_j - X_j - bblockC_j }
                          =>
        { bblockC_i - X_i - X_j - bblockC_j }
        '''
        # Empty container checks
        if bblockC_j.n_particles == 0:             # If struc2 is empty (has no particles)
            return copy.deepcopy(bblockC_i)                 #   simply return unchanged current container
        if bblockC_i.n_particles == 0:              # If struc1 (this struc) is empty (has no particles)
            return copy.deepcopy(bblockC_j)                 #   simply return unchanged current container

        # If no tag for the new structure is provided c
        # ombine the tag of i and j
        if( tag == "blank" ):
            tag = bblockC_i.tag +  bblockC_j.tag
        
        Xn_i = 0  # Number of term atom in neighbor list of cap atom 
        Xn_j = 0  # Number of term atom in neighbor list of cap atom
        #
        # Make copies of containers to modify 
        bbC_i = copy.deepcopy(bblockC_i)
        bbC_i.tag = tag
        bbC_j = copy.deepcopy(bblockC_j)
        #
        # Record attachment 
        attachment_i = Attachment(tag,bblockC_i.tag,bblockC_j.tag,bbid_i,n_i,bbid_j,n_j)
        bbC_i.attachments.append(attachment_i)
        #
        # Find keys of attachment points 
        Rkey_i,Xkey_i = bbC_i.get_rsite(bbid_i,n_i,Xn_i)
        Rkey_j,Xkey_j = bbC_j.get_rsite(bbid_j,n_j,Xn_j)
        #
        # Sum charges of particles to be removed into attachment points
        bbC_i.sum_charge(Xkey_i,Rkey_i)
        bbC_j.sum_charge(Xkey_j,Rkey_j)
        # Align building blocks along bonds of attachment atoms 
        #
        bbC_i.align_bond(Rkey_i,Xkey_i)
        bbC_i.shift_pos(-1.0*bbC_i.positions[Xkey_i] )
        bbC_j.align_bond(Xkey_j,Rkey_j)
        #
        # Shift  building block j to correct bond length 
        radii_i = bbC_i.particles[Xkey_i].bonded_radius
        radii_j = bbC_j.particles[Xkey_j].bonded_radius
        bond_vec = np.array([radii_i + radii_j,0.0,0.0])
        bbC_j.shift_pos(-1.0*bond_vec)
        #
        # Remove atoms at R in building block i
        bbC_i.del_particle(Rkey_i)
        Xo_i = bbC_i.keyupdate[Xkey_i]
        bbC_j.del_particle(Rkey_j)
        Xo_j = bbC_j.keyupdate[Xkey_j]
        #
        # Add j to i
        bbC_i += bbC_j
        #
        # Get updated atom keys to form bond 
        Xp_i = Xo_i                   
        Xp_j = bbC_i.keyupdate[Xo_j]
        #
        # Create bond between  X_i - X_j
        bond_ij = Bond(Xp_i,Xp_j)
        bbC_i.add_bond(bond_ij)
        #
        # Remake neighbor list based on updated bonds 
        bbC_i.bonded_nblist = NBlist()
        bbC_i.bonded_nblist = bbC_i.guess_nblist(0,radii_buffer=1.25)
        # Update number of attachment points
        bbC_i.find_rsites()
        # Set points of attachments
        bbC_i.Xp_i = Xp_i
        bbC_i.Xp_j = Xp_j
        #
        # Return final structure
        # 
        return bbC_i

