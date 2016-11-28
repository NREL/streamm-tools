"""
This module defines the classes relating to a molecular building block 
"""

__author__ = "Travis W. Kemper"
__version__ = "0.3"
__email__ = "travis.kemper.w@gmail.com"
__status__ = "Beta"



import numpy as np 
import copy, sys

import structure, periodictable, units 

from structure import Atom 
from structure import Bond 
from structure import Container as  structure_Container

import logging
logger = logging.getLogger(__name__)


class Attachment():
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

        
class BBatom(Atom):
    """
    A derived type of particles for atoms in a building block 
    """

    def __init__(self,symbol="X", type="bbatom"):
        structure.Atom.__init__(self,symbol=symbol, type=type)
        self.type = type
        self.properties["bbid"] = ""
        self.properties["cplytag"] =  ""  

    def __str__(self):
        """
        'Magic' method for printng contents of container
        """
        return " %s "%(self.properties["symbol"])

    
class Container(structure_Container):
    """
    Data structure for describing a section of a molecule
       that can be joined together to make a new molecule 
    """    
    def __init__(self,tag="blank", verbose=False):    
        """
        Constructor for a composite structure. 
        """
        structure.Container.__init__(self)
        #super(structure.Container,self).__init__()
        self.tag = tag 
        
        self.n_term = int(0)
        self.n_func = int(0)
        self.n_sub = int(0)
        self.terms = []
        self.funcs = []
        self.subs = []        
        self.attachments = []

        self.properties['bblist'] = ""
        
    def __del__(self):
        del self.n_term
        del self.n_func
        del self.n_sub
        del self.terms
        del self.funcs
        del self.subs
        del self.attachments

    def convert_gaussian(self,calc_gaussian):
        
        '''
        Convert Gaussian simulation to Buildingblock

        Input:
           calc_gaussian (Gaussian) simulation object

        Units:
            Bohr -> angstroms

        '''
        import gaussian
        
        if isinstance(calc_gaussian, gaussian.Gaussian):
            #
            # Convert Structure 
            #
            if( calc_gaussian.strucC.n_particles > 0):            
                # Add lattice to GROMACS object
                matrix_o = calc_gaussian.strucC.lat._matrix 
                matrix_i = self.lat._matrix 
                for m in range(calc_gaussian.strucC.lat.n_dim):
                    for n in range(calc_gaussian.strucC.lat.n_dim):
                        matrix_i[m][n] = units.convert_bohr_ang(matrix_o[m][n] )
                self.lat.set_matrix(matrix_i)
                #
                # Add particles and positions to LAMMPS object 
                #
                for pkey_o, particle_o  in calc_gaussian.strucC.particles.iteritems():
                    particle_i = copy.deepcopy(particle_o)
                    pos_o = calc_gaussian.strucC.positions[pkey_o]       
                    pos_i = [ units.convert_bohr_ang(v_i) for v_i in pos_o ] 
                    self.add_partpos(particle_i,pos_i)

                if( len(calc_gaussian.strucC.bonded_nblist.list) > 0 ):
                    self.bonded_nblist = copy.deepcopy(calc_gaussian.strucC.bonded_nblist)
                if( len(calc_gaussian.strucC.nonbonded_nblist.list) > 0 ):
                    self.nonbonded_nblist = copy.deepcopy(calc_gaussian.strucC.nonbonded_nblist)
                    
                for bkey_o, bond_o  in calc_gaussian.strucC.bonds.iteritems():
                    bond_i = copy.deepcopy(bond_o)
                    self.add_bond(bond_i)
                for akey_o, angle_o in calc_gaussian.strucC.angles.iteritems():
                    angle_i = copy.deepcopy(angle_o)
                    self.add_angle(angle_i)
                for dkey_o,dih_o in calc_gaussian.strucC.dihedrals.iteritems():
                    dih_i = copy.deepcopy(dih_o)
                    self.add_dihedral(dih_i)
                for ikey_o,imp_o in calc_gaussian.strucC.impropers.iteritems():
                    imp_i = copy.deepcopy(imp_o)
                    self.add_improper(imp_i)

                    
    def calc_attachments(self):
        '''
        Count number terminal and functional attachment points
        '''
        self.n_term = 0 
        self.n_func = 0
        
        for pkey_i, particle_i  in self.particles.iteritems():
            if( particle_i.properties["bbid"] == "T" ):
                self.n_term += 1 
            if( particle_i.properties["bbid"] == "R" ):
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
           
            if( particle_i.properties["cplytag"][:8] == "termcap_" ):
                particle_i.properties["bbid"] = "T"
                self.terms.append(pkey_i)
                self.n_term += 1 
            elif( particle_i.properties["cplytag"][:8] == "funccap_" ):
                particle_i.properties["bbid"] = "R"
                self.funcs.append(pkey_i)
                self.n_func += 1 
            elif( particle_i.properties["cplytag"][:6] == "rgcap_" ):
                particle_i.properties["bbid"] = "R"
                self.funcs.append(pkey_i)
                self.n_func += 1
            elif( particle_i.properties["cplytag"][:5] == "fcap_" ):
                particle_i.properties["bbid"] = "S"
                self.subs.append(pkey_i)
                self.n_sub += 1
            elif( particle_i.properties["cplytag"][:5] == "term_" ):
                particle_i.properties["bbid"] = "X"
            elif( particle_i.properties["cplytag"][:5] == "func_" ):
                particle_i.properties["bbid"] = "X"
            else: 
                particle_i.properties["bbid"] = particle_i.properties["cplytag"]
                
 
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
           
            if( particle_i.properties["bbid"] == "T" ):                
                self.properties['term_list'].append(pkey_i)
            elif( particle_i.properties["bbid"] == "R"  ):
                self.properties['func_list'].append(pkey_i)
            elif( particle_i.properties["bbid"] == "S"):
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
            
            if( particle_i.properties["bbid"] == bbid_i ):
                # print " debug self.bonded_nblist.getnbs(pkey_i) ",cnt_i,self.bonded_nblist.getnbs(pkey_i)[Xn_i]
                if( cnt_i == n_i ):
                    Rkey_i = pkey_i
                    Xkey_i = self.bonded_nblist.getnbs(pkey_i)[Xn_i]
                    # R - X - StrucC 
                    return Rkey_i,Xkey_i
                cnt_i += 1

        error_line = " In streamm.buildingblock.find_XR \n"
        error_line += " key n=%d with bbid == %s and neighbor %d in building block %s not found "%(n_i,bbid_i,Xn_i,self.tag)
        sys.exit(error_line)


    def align_yaxis(self, key_k,yangle,debug = False):
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

        if( debug ):
            print " pos_k ",pos_k
            print "     >align_yaxis  pos_k", pos_k
        
        # Get the zy component of position k
        pos_yz = np.array([0.0,pos_k[1],pos_k[2]])
        # Unit vector in the y direction 
        unit_y = np.array( [0.0,1.0,0.0])

        # Rotate structure around x
        # according to angle between position k and y-axis
        # Subtract desired final  angle with the y-axis
        if( debug ):
            print "     >align_yaxis  yangle", yangle
        cos_yz = calc_cos(unit_y,pos_yz)
        angle_yz = np.arccos(cos_yz)
        if( debug ):
            print "     >align_yaxis  angle_yz", angle_yz
        angle_yz += -1.0*yangle
        
        if( debug ):
            print "     >align_yaxis  angle_yz - yangle ", angle_yz
        
        if( pos_yz[2] > 0.0 ):
            # if in the 1st or 2nd quadrents rotate  clockwise
            direction="clockwise"
        elif( pos_yz[2] <= 0.0 ):
            # if in the 3rd or 4th quadrents rotate counter clockwise 
            # of if 180 rotation
            direction="counterclockwise"
        if( debug ):
            print "     >align_yaxis  direction", direction
        self.rotate_yz(angle_yz,direction=direction,verbose=False)
        # if  pos_xy[1] == 0.0 no rotation is needed
        if( debug ):
            pos_k = np.array(self.positions[key_k]  )
            print "     >align_yaxis  pos_k", pos_k
        
        return

    def align_bond(self, key_i, key_j,debug = False):
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

        if( debug ):
            print ">align_bond pos_i pos_j ",pos_i,pos_j

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
        if( debug ):
            print "     >align_bond pos_xz",pos_xz
        
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

        if( debug ):
            print "     >align_bond  after xz rot self.positions[key_j] ",self.positions[key_j] 

        return
    
    def set_dih(self,key_l,key_i,key_j,key_k,angle):
        '''
        Set the dihedral angle between key_l,key_i,key_j,key_k  to angle

        Assume that the bond between key_i,key_j is along the x-axis using
        align_bond(self, key_i, key_j)
        
        '''
        cos_kijl = self.calc_cosdih(key_k,key_i,key_j,key_l)
        
        
    def read_cply(self,cply_file='',debug = False):
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
            tag_s2 = tag_s1[-1]            # take last string 
            self.tag = tag_s2[:-5]         # remove .cply
            self.properties['common_tag'] = tag_s2[:-5]         # remove .cply

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
                        BBatom_i.properties["cplytag"] = str(col[13])

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
                        BBatom_i.properties["cplytag"] = str(col[7])
                        
                    elif (len(col) == 7 ):
                        BBatom_i.properties["charge"] = float(col[4])
                        BBatom_i.properties["residue"] = int(col[5])
                        BBatom_i.properties["resname"] = str(col[6])

                    elif (len(col) == 5 ):
                        BBatom_i.properties["cplytag"] = str(col[4])

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
                        self.properties['bblist'] = str(col[3:])

                # Key word search 
                elif( len(col) > 0 ):
                    if( str(col[0]) == 'lattice' ):
                        read_lattice = True
                        lv_cnt = 0                        
        if( debug ):
            sys.exit("debug in read_cply is True ")
        
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
        cply_line += "type %s %s %s \n"%(self.properties['moltype'],self.properties['backbone'],self.properties['bblist'])
        
        if(write_ff ):
            cply_line += "# atomic_symb ,float(pos_i[0]), float(pos_i[1]),float(pos_i[2]),label,fftype,ptclObj.mass,charge,qgroup,ring,residue,resname, mol, cplytag \n"
        else:
            cply_line += "# atomic_symb ,float(pos_i[0]), float(pos_i[1]),float(pos_i[2]), cplytag \n"
            
        F.write(cply_line)
        for pkey_i, particle_i  in self.particles.iteritems():
            pos_i = self.positions[pkey_i]
            atomic_symb = particle_i.properties['symbol']
            bbid = particle_i.properties["bbid"]
            cplytag = particle_i.properties["cplytag"]
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
    
    def prepattach(self,bbid_i="R",n_i=0,Xn_i=0,dir=1,yangle=0.0,debug = False):
        '''
        Preppair a building block to be attached to another

        This is handy for batch attachments to not repeat the intial steps in the attach function

        Arguments:
            yangle (float) angle in radians
            
        '''
    
        if( debug ):                
            print " >prepattach ",bbid_i,n_i,Xn_i
            print "self.positions 9247502 ",bb_prepped.positions
        
        bb_prepped = copy.deepcopy(self)
        if( bb_prepped.n_particles > 1 ):
            
            # Find keys of attachment points
            #   j - i - Struc
            # Rkey_i,Xkey_i
            key_j,key_i = bb_prepped.find_XR(bbid_i,n_i,Xn_i)
            # Align building blocks along bonds of attachment atoms
            if( dir == 1 ):
                bb_prepped.align_bond(key_i,key_j,debug = debug)
                if( debug ):                
                    pos_i = bb_prepped.positions[key_i]
                    print " > prepattach  align_bond pos_i ",dir,pos_i,key_i
            elif( dir == -1 ):
                bb_prepped.align_bond(key_j,key_i,debug = debug)
                # 
                pos_i = bb_prepped.positions[key_i]
                bb_prepped.shift_pos(-1.0*pos_i)
                if( debug ):                
                    print " > prepattach  align_bond pos_i ",dir,pos_i,key_i
            else:
                error_line = " dir is not 1 or -1 "
                error_line += " in streamm.buildingblock.attachprep "
                print error_line

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


def checkprep(bbC_i,bbC_j,covbuffer=1.5,debug=False):
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
    debug = False
    if( debug ):
        print " >checkprep debug on "
        
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
        if( pkey_i != Xo_i and particle_i.properties['number'] != 1 ):
            for pkey_j, particle_j  in bbC_j.particles.iteritems():
                if( pkey_j != Xo_j and particle_j.properties['number'] != 1 ):
                    dij = nd_ij[pkey_i][pkey_j]
                    radii_i = bbC_i.particles[pkey_i].properties["cov_radii"]
                    radii_j = bbC_j.particles[pkey_j].properties["cov_radii"]
                    cut_off = radii_i +radii_j
                    cut_off = cut_off*covbuffer
                    if( dij < cut_off ):
                        if( debug ):
                            log_line = "      >checkprep \n"
                            log_line += "           particle i %s %d \n"%(particle_i.properties['symbol'],pkey_i)
                            log_line += "           particle j %s %d \n"%(particle_j.properties['symbol'],pkey_j)
                            log_line += "           radii_i %f \n"%(radii_i)
                            log_line += "           radii_j %f \n"%(radii_j)
                            log_line += "           dij %f \n"%(dij)
                            print log_line
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
    radii_i = bbC_i.particles[Xo_i].properties["cov_radii"]
    radii_j = bbC_j.particles[Xo_j].properties["cov_radii"]
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
    
def attachprep(bbC_i,bbC_j,debug=False):
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
    bond_ij = structure.Bond(Xo_i,Xp_j)
    bond_ij.properties['bondorder'] = 1
    print ">attachprep set bondorder ",Xo_i,Xp_j,bond_ij.properties['bondorder']
    bbC_i.add_bond(bond_ij)
    #
    # Remake neighbor list based on updated bonds 
    bbC_i.bonded_nblist = structure.NBlist() 
    bbC_i.bonded_nblist.build_nblist(bbC_i.particles,bbC_i.bonds )
    # Update number of attachment points
    bbC_i.calc_attachments()                

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
        Rkey_i,Xkey_i = bbC_i.find_XR(bbid_i,n_i,Xn_i)
        Rkey_j,Xkey_j = bbC_j.find_XR(bbid_j,n_j,Xn_j)
        #
        # Sum charges of particles to be removed into attachment points
        bbC_i.sum_prop("charge",Xkey_i,Rkey_i)
        bbC_j.sum_prop("charge",Xkey_j,Rkey_j)
        # Align building blocks along bonds of attachment atoms 
        #
        bbC_i.align_bond(Rkey_i,Xkey_i)
        bbC_i.shift_pos(-1.0*bbC_i.positions[Xkey_i] )
        bbC_j.align_bond(Xkey_j,Rkey_j)
        #
        # Shift  building block j to correct bond length 
        radii_i = bbC_i.particles[Xkey_i].properties["cov_radii"]
        radii_j = bbC_j.particles[Xkey_j].properties["cov_radii"]
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
        bond_ij = structure.Bond(Xp_i,Xp_j)
        bbC_i.add_bond(bond_ij)
        #
        # Remake neighbor list based on updated bonds 
        bbC_i.bonded_nblist = structure.NBlist() 
        bbC_i.bonded_nblist.build_nblist(bbC_i.particles,bbC_i.bonds )
        # Update number of attachment points
        bbC_i.calc_attachments()
        # Set points of attachments
        bbC_i.Xp_i = Xp_i
        bbC_i.Xp_j = Xp_j
                      
        #  Return final structure 
        return bbC_i


def calc_cos(v_i,v_j):
    """
    Calculate the cos theta between two vectors 
    """
    return np.dot(v_i/np.linalg.norm(v_i),v_j/np.linalg.norm(v_j))
