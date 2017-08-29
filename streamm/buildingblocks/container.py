.bbid# coding: utf-8
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
    
    # Import streamm dependencies
    from streamm.structure.containers import Container as strucCont
    # from streamm.structure.container import Container
    from streamm.structure.lattices import Lattice 
    from streamm.structure.nblists import NBlist 
    
    from streamm.structure.atoms import Atom
    from streamm.structure.bonds import Bond 
    from streamm.structure.angles import Angle
    from streamm.structure.dihedrals import Dihedral
    from streamm.structure.impropers import Improper

except:
    print("streamm is not installed test will use relative path")
    import sys, os
    rel_path = os.path.join(os.path.dirname(__file__),'..','structure')
    print("rel_path {}".format(rel_path))
    sys.path.append(rel_path)
    
    from containers import Container as strucCont
    from particletype import Particletype
    from bondtype import Bondtype 
    from angletype import Angletype
    from dihtype import Dihtype
    from imptype import Imptype


class Attachment(object):
    '''
    Object to recode the attachment of two building block objects
    '''
    def __init__(self,tag_ij,tag_i,tag_j,bbid_i,n_i,bbid_j,n_j):

        self.tag_ij = tag_ij
        self.tag_i = tag_i
        self.tag_j = tag_j
        self.bbid_i = bbid_i
        self.n_i = n_i
        self.bbid_j = bbid_j
        self.n_j = n_j
    
    def __del__(self):
        """
        'Magic' method for deleting contents of container
        """
        del self.tag_ij
        del self.tag_i
        del self.tag_j
        del self.bbid_i
        del self.n_i
        del self.bbid_j
        del self.n_j

    def __str__(self):
        """
        'Magic' method for printng contents of container
        """
        return " %s ( %s %d ) +  %s ( %s %d ) -> %s "%(self.tag_i,self.bbid_i,self.n_i,self.tag_j,self.bbid_j,self.n_j,self.tag_ij)
    
class Container(strucCont):
    """
    Data structure for describing a section of a molecule
       that can be joined together to make a new molecule 
    """    
    def __init__(self,tag="blank", verbose=False):    
        """
        Constructor for a composite structure. 
        """
        strucCont.__init__(self)
        #super(structure.Container,self).__init__()
        self.tag = tag 
        
        self.n_term = int(0)
        self.n_func = int(0)
        self.n_sub = int(0)
        self.terms = []
        self.funcs = []
        self.subs = []        
        self.attachments = []

        self.bblist = ""
        
        
        
    def __del__(self):
        del self.n_term
        del self.n_func
        del self.n_sub
        del self.terms
        del self.funcs
        del self.subs
        del self.attachments
        del self.bblist
        
    def calc_attachments(self):
        '''
        Count number terminal and functional attachment points
        '''
        self.n_term = 0 
        self.n_func = 0
        
        for pkey_i, particle_i  in self.particles.iteritems():
            if( particle_i.bbid == "T" ):
                self.n_term += 1 
            if( particle_i.bbid == "R" ):
                self.n_func += 1 
                
    def parse_cplytag(self):
        '''
        Define bbid pased on <=v0.2 cplytag

        bb - X - T  - terminal site for polymerization 
        bb - X - R  - functionalizable site for functional groups 
        bb - X - S  - substitutable site for atomic substitutions 

        this should be flipped ...
        but it makes it easier to id the H since it's only bonded to 1 thing
         
        '''
        self.n_term = 0 
        self.n_func = 0
        self.n_sub = 0

        self.terms = []
        self.funcs = []
        self.subs = []
                
        for pkey_i, particle_i  in self.particles.iteritems():
           
            if( particle_i.cplytag[:8] == "termcap_" ):
                particle_i.bbid = "T"
                self.terms.append(pkey_i)
                self.n_term += 1 
            elif( particle_i.cplytag[:8] == "funccap_" ):
                particle_i.bbid = "R"
                self.funcs.append(pkey_i)
                self.n_func += 1 
            elif( particle_i.cplytag[:6] == "rgcap_" ):
                particle_i.bbid = "R"
                self.funcs.append(pkey_i)
                self.n_func += 1
            elif( particle_i.cplytag[:5] == "fcap_" ):
                particle_i.bbid = "S"
                self.subs.append(pkey_i)
                self.n_sub += 1
            elif( particle_i.cplytag[:5] == "term_" ):
                particle_i.bbid = "X"
            elif( particle_i.cplytag[:5] == "func_" ):
                particle_i.bbid = "X"
            else: 
                particle_i.bbid = particle_i.cplytag
                
 
    def proc_bbid(self):
        '''
        Define bbid pased on <=v0.2 cplytag

        bb - X - T  - terminal site for polymerization 
        bb - X - R  - functionalizable site for functional groups 
        bb - X - S  - substitutable site for atomic substitutions 

        this should be flipped ...
        but it makes it easier to id the H since it's only bonded to 1 thing
         
        '''

        self.properties['term_list'] = []
        self.properties['func_list']  = []
        self.properties['sub_list']  = []
            
        for pkey_i, particle_i  in self.particles.iteritems():
           
            if( particle_i.bbid == "T" ):                
                self.properties['term_list'].append(pkey_i)
            elif( particle_i.bbid == "R"  ):
                self.properties['func_list'].append(pkey_i)
            elif( particle_i.bbid == "S"):
                self.properties['sub_list'].append(pkey_i)
                
        return 
        

    def find_XR(self,bbid_i,n_i,Xn_i):
        '''
        Find keys of terminal atom X and it's cap R
        
        Arguments:
            bbid_i   (str) Connecting atom bbid in container 1 
            n_i      (int) number of connection to used in  container 1 
            Xn_i     (int) nieghbor number of R atom

        Return:
            Rkey  - the particle to be replaced
            XKey  - the nieghbor Xn_i of the particle to be replaced 
        '''
        cnt_i = 0

        # print " find_XR self.tag ",self.tag
        
        for pkey_i, particle_i  in self.particles.iteritems():
            
            if( particle_i.bbid == bbid_i ):
                # print " debug self.bonded_nblist.getnbs(pkey_i) ",cnt_i,self.bonded_nblist.getnbs(pkey_i)[Xn_i]
                if( cnt_i == n_i ):
                    Rkey_i = pkey_i
                    Xkey_i = self.bonded_nblist.getnbs(pkey_i)[Xn_i]
                    # R - X - StrucC 
                    return Rkey_i,Xkey_i
                cnt_i += 1
        error_line = "Particle key n=%d with bbid = %s and neighbor %d in building block %s not found "%(n_i,bbid_i,Xn_i,self.tag)
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
    
    def set_dih(self,key_l,key_i,key_j,key_k,angle):
        '''
        Set the dihedral angle between key_l,key_i,key_j,key_k  to angle

        Assume that the bond between key_i,key_j is along the x-axis using
        align_bond(self, key_i, key_j)
        
        '''
        cos_kijl = self.calc_cosdih(key_k,key_i,key_j,key_l)
        
        
    def read_cply(self,cply_file=''):
        """
        Read cply file
        """
        if( len(cply_file) == 0 ):
            cply_file = "%s.cply"%(self.tag)
        read_lattice = False
        read_attachments = False

        # Set tag of structure based on file
        if( self.tag == "blank" ):
            tag_s1 = cply_file.split("/")   # split string based on dir 
            tag_s2 = tag_s1[-1]             # take last string 
            self.tag = tag_s2[:-5]          # remove .cply
            self.properties['common_tag'] = tag_s2[:-5]         # remove .cply
        #
        with open(cply_file) as f:
            for line in f:
                col = line.split()
                if( read_lattice ):
                    matrix = self.lat._matrix
                    matrix[lv_cnt][0] = float( col[0] )
                    matrix[lv_cnt][1] = float( col[1] )
                    matrix[lv_cnt][2] = float( col[2] )
                    if( lv_cnt == 2 ):
                        self.lat.set_matrix(matrix)
                        read_lattice = False 
                    lv_cnt += 1
                    
                elif( len(col) >= 4  and col[0] != "name" and col[0] != "bond"  and col[0] != "attachment"  and col[0] != "replication"  and col[0] != "#" and col[0] != "type"  and col[0] != "func" ):

                    BBatom_i = BBatom(str(col[0]))
                    
                    pos_i = [ float(col[1]),float(col[2]),float(col[3]) ]

                    if (len(col) >= 14 ):
                        BBatom_i.properties["label"] = str(col[4])
                        BBatom_i.properties["fftype"] = str(col[5])
                        BBatom_i.properties["mass"] = float(col[6])                        
                        BBatom_i.properties["charge"]  = float(col[7])
                        BBatom_i.properties["qgroup"] = int(col[8])                        
                        BBatom_i.properties["ring"] = int(col[9])                        
                        BBatom_i.properties["residue"] = int(col[10])
                        BBatom_i.properties["resname"] = str(col[11])                       
                        BBatom_i.properties["mol"] = int(col[12])                        
                        BBatom_i.cplytag = str(col[13])

                    elif (len(col) == 13 ):
                        BBatom_i.properties["label"] = str(col[4])
                        BBatom_i.properties["fftype"] = str(col[5])
                        BBatom_i.properties["mass"] = float(col[6])                        
                        BBatom_i.properties["charge"] = float(col[7])
                        BBatom_i.properties["qgroup"] = int(col[8])                        
                        BBatom_i.properties["ring"] = int(col[9])                        
                        BBatom_i.properties["residue"] = int(col[10])
                        BBatom_i.properties["resname"] = str(col[11])                       
                        BBatom_i.properties["mol"] = int(col[12])

                    elif (len(col) == 8 ):
                        BBatom_i.properties["residue"] = int(col[5])
                        BBatom_i.properties["resname"] = str(col[6])                        
                        BBatom_i.cplytag = str(col[7])
                        
                    elif (len(col) == 7 ):
                        BBatom_i.properties["charge"] = float(col[4])
                        BBatom_i.properties["residue"] = int(col[5])
                        BBatom_i.properties["resname"] = str(col[6])

                    elif (len(col) == 5 ):
                        BBatom_i.cplytag = str(col[4])

                    self.add_partpos(BBatom_i,pos_i, deepcopy = True)
                   
                elif(len(col) >= 3 ):
                    if( col[0] == "bond"):
                        pkey1 = int(col[1]) - 1
                        pkey2 = int(col[2]) - 1
                        bond_i = Bond(pkey1,pkey2)
                        if( len(col) >= 4 ):
                            bond_i.border =  int(col[3]) 
                        #print "process_line bond line ",col
                        self.add_bond(bond_i)
                    if( col[0] == "name"):
                        self.tag = str(col[1])
                        self.properties['common_tag'] = str(col[1])
                        self.properties['deptag'] = str(col[2])
                        self.properties['ctag'] = str(col[2])
                        self.properties['name'] = str(col[3])
                        self.properties['IUPAC'] = ''
                        for s in col[4:]:
                            self.properties['IUPAC'] += str(s)
                    if( col[0] == "type"):
                        self.properties['moltype'] = str(col[1])
                        self.properties['backbone'] = str(col[2])
                        self.bblist = str(col[3:])

                # Key word search 
                elif( len(col) > 0 ):
                    if( str(col[0]) == 'lattice' ):
                        read_lattice = True
                        lv_cnt = 0                        
          
        # Convert old tags to new bbid 
        self.parse_cplytag()
        # Build neighbor list of bonded atoms 
        if( self.n_bonds == 0 ):
            # If no bonds guess based on radii 
            self.bonded_nblist.guess_nblist(self.lat,self.particles,self.positions,"cov_radii",radii_buffer=1.25)
            # Build bonds from nblist for reference 
            self.bonded_bonds()
        else:
            self.bonded_nblist.build_nblist(self.particles,self.bonds )

    def write_cply(self,cply_file='', write_ff=True,write_bonds=True,write_latvec=True,write_attachments=True):
        """
        Write cply file
        """
        if( len(cply_file) == 0 ):
            cply_file = "%s.cply"%(self.tag)

        F = open(cply_file,'w')
        cply_line = "name %s %s %s %s \n"%(self.tag,self.properties['deptag'],self.properties['name'],self.properties['IUPAC'])
        cply_line += "type %s %s %s \n"%(self.properties['moltype'],self.properties['backbone'],self.bblist)
        
        if(write_ff ):
            cply_line += "# atomic_symb ,float(pos_i[0]), float(pos_i[1]),float(pos_i[2]),label,fftype,ptclObj.mass,charge,qgroup,ring,residue,resname, mol, cplytag \n"
        else:
            cply_line += "# atomic_symb ,float(pos_i[0]), float(pos_i[1]),float(pos_i[2]), cplytag \n"
            
        F.write(cply_line)
        for pkey_i, particle_i  in self.particles.iteritems():
            pos_i = self.positions[pkey_i]
            atomic_symb = particle_i.properties['symbol']
            bbid = particle_i.bbid
            cplytag = particle_i.cplytag
            if(write_ff ):
                mass = particle_i.properties["mass"]
                charge = particle_i.properties["charge"]
                residue = particle_i.properties["residue"]
                resname = particle_i.properties["resname"]
                label = particle_i.properties["label"]
                fftype = particle_i.properties["fftype"]
                qgroup = particle_i.properties["qgroup"]
                mol = particle_i.properties["mol"]
                ring = particle_i.properties["ring"]
                cply_line =  " %5s %16.8f %16.8f %16.8f %s %s %12.8f %12.8f  %d %d %d %s  %d %s \n"%(atomic_symb ,float(pos_i[0]), float(pos_i[1]),float(pos_i[2]),label,fftype,mass,charge,qgroup,ring,residue,resname, mol, bbid )
                #cply_line =  "%d %5s %16.8f %16.8f %16.8f %s %s %12.8f %12.8f  %d %d %d %s  %d %s \n"%(pkey_i,atomic_symb ,float(pos_i[0]), float(pos_i[1]),float(pos_i[2]),label,fftype,mass,charge,qgroup,ring,residue,resname, mol, cplytag )
            else:
                cply_line =  " %5s %16.8f %16.8f %16.8f %s \n"%(atomic_symb ,float(pos_i[0]), float(pos_i[1]),float(pos_i[2]), cplytag )
            F.write(cply_line)

        if( self.lat._matrix[0][0] == 0.0 ):
            write_latvec = False 
        if( write_latvec ):
            F.write(" lattice \n")
            for d in range(3):
                cply_line = " %16.8f %16.8f %16.8f \n"%( self.lat._matrix[d][0],self.lat._matrix[d][1],self.lat._matrix[d][2])
                F.write(cply_line)
            

        if( write_bonds ):
            if( self.n_bonds > 0 ):
                for bkey_i, bond_i  in self.bonds.iteritems():
                    b_i = bond_i.pkey1 + 1 
                    b_j = bond_i.pkey2 + 1
                    # F.write("  bond %d %d  %d \n"%(b_i,b_j,bond_i.properties['order']))
                    F.write("  bond %d %d  \n"%(b_i,b_j))

        if( write_attachments ):
            for attachment_i in self.attachments :
                cply_line = " attachment "
                cply_line += str(attachment_i) +"\n"
                F.write(cply_line)
            for replication_i in self.replications :
                cply_line = " replication "
                cply_line += str(replication_i) +"\n"
                F.write(cply_line)                
        F.close()

        return

    def getSubStructure(self,pkeys,tag="blank"):
        """
        This is redundant to getSubStructure in structure.py
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
            new_strucC.bonded_nblist = structure.NBlist() 
            for pkey_i in pkeys:
                new_strucC.bonded_nblist.index.append(new_strucC.bonded_nblist.cnt + 1)
                for pkey_j in self.bonded_nblist.getnbs(pkey_i):
                    if( pkey_j in pkeys ):
                        new_strucC.bonded_nblist.cnt += 1 
                        new_strucC.bonded_nblist.list.append(key_update[pkey_j])

            new_strucC.bonded_nblist.index.append(new_strucC.bonded_nblist.cnt + 1)
            new_strucC.bonded_bonds()

        return new_strucC
    
    def prepattach(self,bbid_i="R",n_i=0,Xn_i=0,dir=1,yangle=0.0):
        '''
        Preppair a building block to be attached to another

        This is handy for batch attachments to not repeat the intial steps in the attach function

        Arguments:
            yangle (float) angle in radians
            
        '''
    
        logger.debug(".prepattach {} {} {}".format(bbid_i,n_i,Xn_i))
        
        if( dir != 1 and dir != -1 ):
            error_line = " dir is not 1 or -1 "
            raise  ValueError(error_line)
            
        
        bb_prepped = copy.deepcopy(self)
        if( bb_prepped.n_particles > 1 ):
            
            # Find keys of attachment points
            #   j - i - Struc
            # Rkey_i,Xkey_i
            key_j,key_i = bb_prepped.find_XR(bbid_i,n_i,Xn_i)
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
                    if( particle_k.properties['number'] != 1 ):                        
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

