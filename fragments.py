"""
This module defines classes and functions for manipulating molecular "fragments".
"""

import os,sys,fileinput,random
# from math import *
import math

from basics import *

class TransRot:
    """
    Data for storing translation and rotation
    """
    def __init__(self):
        self.trans = None
        self.rot = None
    def __init__(self, T,R):
        self.trans = T
        self.rot = R


def get_axis_angle_rotation(n, theta):
    """
    return 3x3 rotation matrix for rotation of theta around axis n (vector)
    Key? for rotating fragments
    Mostly the self.rot matrix in class above is generated here
    """
    # c = cos(theta)
    # s = sin(theta)
    c = math.cos(theta)
    s = math.sin(theta)
    t = 1.0 - c

    m00 = c + n[0]*n[0]*t
    m11 = c + n[1]*n[1]*t
    m22 = c + n[2]*n[2]*t

    tmp1 = n[0]*n[1]*t
    tmp2 = n[2]*s
    m10 = tmp1 + tmp2
    m01 = tmp1 - tmp2
    tmp1 = n[0]*n[2]*t
    tmp2 = n[1]*s
    m20 = tmp1 - tmp2
    m02 = tmp1 + tmp2
    tmp1 = n[1]*n[2]*t
    tmp2 = n[0]*s
    m21 = tmp1 + tmp2
    m12 = tmp1 - tmp2

    R = [[m00, m01, m02],[m10,m11,m12],[m20,m21,m22]]
    return R

def get_2pt_transform(x1,x2,y1,y2):
    """

    Args:
        x1 (list pts) origin for x-vector
        x2 () other for x-vector
        y1 () origin for y-vector
    Return:
        TransRot object
    """
    trans = vec_sub(x1,y1)  # vec calls in basics
    z = vec_sub(x2,x1)
    w = vec_sub(y2,y1)
    # theta = acos(vec_dot(z,w)/(vec_norm(z) * vec_norm(w)))
    theta = math.acos(vec_dot(z,w)/(vec_norm(z) * vec_norm(w)))
    n = vec_cross(w,z)
    n = vec_normalize(n)
    #    print "n, theta", n, theta
    R = get_axis_angle_rotation(n, theta)
    return TransRot(trans, R)


def apply_transforms(TR, R2, a, pos0):
    """
    Only called in fragments.py
    Args:
        TR: 
        R2: 
        a
        pos0
    """
    # apply shift
    tmp = vec_add(a.pos, TR.trans)
    # apply axis angle rotation 1, with pos0 as the origin
    tmp = vec_sub(tmp, pos0)
    tmp = mat_vec_mult(TR.rot, tmp)
    tmp = vec_add(tmp, pos0)
    # apply axis angle rotation 2, with pos0 as the origin
    tmp = vec_sub(tmp, pos0)
    tmp = mat_vec_mult(R2, tmp)
    ta = vec_add(tmp, pos0)

    return ta


def get_normal_from_PCA(atompos):
    """
    Get approx normal to points(atoms) in structure

    first we form A = X^tX where X is matrix whose rows are positions in 3D
    then we get its eigenvalues
    The normal is the eigenvector for the smallest eigenvalue.
    Get this from Cayley-Hamilton

    Args:
        atompos: list of atoms positions [ [x1,y1,z1], [x2,y2,z2].....]
    Returns: 
        Approx normal vector
    """

    # Given an real symmetric 3x3 matrix A, compute the eigenvalues
    from math import sqrt, acos, cos, pi

    A = matT_mat_mult(atompos, atompos)
    
    p = A[0][1]**2 + A[0][2]**2 + A[1][2]**2
    if (p == 0): 
        # A is diagonal.
        eig1 = A[0][0]
        eig2 = A[1][1]
        eig3 = A[2][2]
    else:
        q = mat3_trace(A)/3.
        p = (A[0][0] - q)**2 + (A[1][1] - q)**2 + (A[2][2] - q)**2 + 2 * p
        p = sqrt(p / 6.)
        B = [[1/p * (A[i][j] - q *(i==j)) for i in range(3)] for j in range(3)]
        # print B
        # B = (1 / p) * (A - q * I)       # I is the identity matrix
        r = mat3_det(B) / 2.
 
        # In exact arithmetic for a symmetric matrix  -1 <= r <= 1
        # but computation error can leave it slightly outside this range.
        if (r <= -1): 
            phi = pi / 3.
        elif (r >= 1):
            phi = 0
        else:
            phi = acos(r) / 3.
 
        # print r,phi
        # the eigenvalues satisfy eig3 <= eig2 <= eig1
        eig1 = q + 2 * p * cos(phi)
        eig3 = q + 2 * p * cos(phi + pi * (2/3.))
        eig2 = 3 * q - eig1 - eig3     # since trace(A) = eig1 + eig2 + eig3

#   print eig1, eig2, eig3
    # Now use these to get normal (eig3 is the smallest, that corresponds to the normal)
    # M = (A-eig1*I)(A-eig2*I)
    M1 = [[A[i][j] - eig1*(i==j) for i in range(3)] for j in range(3)]
    M2 = [[A[i][j] - eig2*(i==j) for i in range(3)] for j in range(3)]
    M = mat_mat_mult(M1,M2)
    foundv = False
    j = 0;
    while (not foundv):
        v = [M[i][j] for i in range(3)]
        if (vec_norm(v) > 1e-10):
            foundv = True
        else:
            j+=1
            assert(j<3),"nonzero normal principal component not found"
    v = vec_normalize(v)
#    print M
#    print v
    return v


class Link:
    """
    type can be term,func,rg,termcap,funccap,rgcap
    These are tags in the cply files
    This object is created for each 'tag (3rd column) in cply file
    Parsed cply file will have atom positions and these Link objects (?)

    mypos (redundant, should delete TODO) is index of the atom this link is associated with.
    targetpos is index of atom this link is, eg., a cap for.  Keeping track of which atoms are caps
    for which atoms is the main job of this class.  The cap atoms (e.g. H) get deleted when we join structures,
    main use of the Link class is seen in function get_connectors().  Other references are just keeping the links
    up to date.
    In a real language these would be pointers, but that's a different headache.
    """

    def __init__(self):
        self.type = None
        self.mypos = -1       # Atom index
        self.targetpos = -1   # Atom index

    def __init__(self, type, a1, a2):
        self.type = type       # term, func, rg etc
        self.mypos = a1
        self.targetpos = a2

    def shift(self, n):
        """
        Called when atom is deleted with a lower index
        """
        self.mypos += n
        self.targetpos += n

    def dump(self, return_str=False):
        s = "link object: type=%s mypos=%d targetpos=%d"  % (self.type, self.mypos, self.targetpos)
        if (return_str):
            return s
        else:
            print s
    

class Atom:
    """ 
    The atom class keeps travk of its position, element, and
    who its neighbors are.  Also may have a "link" object which describes the identification of a part of 
    one structure with that of another
    """
    def __init__(self):
        self.pos = [0,0,0]
        self.elmnt = "NA"
        self.nb1 = []
        self.nb2 = []
        self.nb = []
        self.nball = []
        self.link = None

        self.ctype = []
        self.resnumb = []
        self.restype = []
        self.q = []


    def __init__(self, elmnt, pos):
        self.pos = pos
        self.elmnt = elmnt
        self.nb1 = []
        self.nb2 = []
        self.nb = []
        self.nball = []
        self.link = None
        
        # make_ff TK
        self.id = ""
        self.ctype = "NA"
        self.resnumb =  1
        self.restype = "NA"
        self.q = 0.0

    def dist(self, other):
        """
        Cartesian dist (local convien)
        """
        return vec_norm(vec_sub(self.pos, other.pos))






def get_basis_str_from_str(accuracy):
    if (accuracy == "high"):
        return "6-31+g(d,p)"
    elif (accuracy == "low"):
        return "6-31g(d)"
    else:
        return "6-31g(d)"
        

def get_basis_str(accuracy):
    """accuracy=0 implies a 3-21g basis
    accuracy=1 implies a 6-31g(d) basis
    accuracy=2 implies a 6-31+g(d) basis
    accuracy=3 implies a 6-31+g(d,p) basis
    """

    s = "6-31g(d)"
    if (accuracy==0):
        s = "3-21g"
    elif (accuracy==1):
        s = "6-31g(d)"
    elif (accuracy==2): 
        s = "6-31+g(d)"
    elif (accuracy==3): 
        s = "6-31+g(d,p)"
    else:
        print "WARNING: unknown accuracy value %d, using default basis %s" % (accuracy, s)

    return s
        

class Structure:
    """
    Base class 
    keeps track of a list of atoms (without connectivity)
    Container for Atom objects
    
    Can calculate intra-building block connectivity through distance
    """

    def __init__(self, struct=None):
        import copy
        self.atoms = []
        if (struct != None):
            self.atoms = [copy.deepcopy(struct.atoms[i]) for i in range(0,len(struct.atoms))]
            self.find_neighbors()
            if (struct.fname != None):
                self.fname = struct.fname

    def write_xyz(self, fname, tag="A Structure"):
        """
        IO writes out .xyz file
        """
        # write data to xyz file for VMD
        f = file(fname, "w")
        f.write("%d\n%s\n" % (len(self.atoms), tag))
        for a in self.atoms:
            f.write("%s %f %f %f\n" % (a.elmnt, a.pos[0], a.pos[1], a.pos[2]))
        f.close()

    def write_com(self, com_template_name, xyz_name, job_name, basis, nstates):
        """
        Writes out complete Gaussian input file (structure part)
        """
        from string import replace
        f = file(com_template_name)
        com_templ = f.read()
        f.close()
        f = file(xyz_name)
        atomlist = f.readlines()[2:]
        atoms = ""
        for at in atomlist:
            atoms += at
        #    print "*********ATOMS********"
        #    print atoms
        #    print  "***************"
        f.close()
        com_templ = replace(com_templ, "<structure_name>", job_name)
        com_templ = replace(com_templ, "<atoms>", atoms)
        com_templ = replace(com_templ, "<basis>", basis)
        com_templ = replace(com_templ, "<nstates>", "%d" % nstates)
        com_name = xyz_name
        com_name = replace(com_name, "xyz", "com")
        f = file(com_name, "w")
        f.write(com_templ)
        f.close()


    def write_com_restart(self, com_template_name, xyz_name, job_name, basis, nstates):
        from string import replace

        # Parse out extension
        tmp=com_template_name.split('.')
        tmp2=tmp[-1]                     # picks out eg 'r1'
        extensionName=tmp2 + ".com"      # constructs eg "r1.com"

        f = file(com_template_name)
        com_templ = f.read()
        f.close()
        com_templ = replace(com_templ, "<structure_name>", job_name)
        com_templ = replace(com_templ, "<basis>", basis)
        com_templ = replace(com_templ, "<nstates>", "%d" % nstates)
        com_name = xyz_name
        com_name = replace(com_name, "xyz", extensionName)
        f = file(com_name, "w")
        f.write(com_templ)
        f.close()


    def write_slurm(self, slurm_template_name, xyz_name, job_name):
        from string import replace
        f = file(slurm_template_name)
        templ = f.read()
        f.close()
        templ = replace(templ, "<job_name>", job_name)
        slurm_name = replace(xyz_name, "xyz", "slurm")
        f = file(slurm_name, "w")
        f.write(templ)
        f.close()
    
    def write_pbs(self, pbs_template_name, xyz_name, job_name):
        from string import replace
        f = file(pbs_template_name)
        templ = f.read()
        f.close()
        templ = replace(templ, "<job_name>", job_name)
        pbs_name = replace(xyz_name, "xyz", "pbs")
        f = file(pbs_name, "w")
        f.write(templ)
        f.close()
    
    def write_meta(self, d,  xyz_name):
        from string import replace
        import json
        meta_name = replace(xyz_name, "xyz", "meta")
        json_name = replace(xyz_name, "xyz", "json")
        f = file(meta_name, "w")
        for key,value in d.items():
            f.write("%s = %s\n" % (key,value))
        f.close()

        f = open(json_name, 'w')
        json.dump(d, f, indent=2)
        f.close()


    def read(self, fname):
        """
        Reads cply file and writes xyz as debug
        """

        self.fname = fname
        # read from piece of Gaussian input file (lines are like xyz: "atom x y z")
        cnt = 0
        for ln in fileinput.input(fname):
            self.process_line(ln, cnt)
            cnt += 1
        # self.dump()
        self.write_xyz("%s.xyz" % self.fname)


    def process_line(self, ln, cnt):
        """
        Parses one line of cply files
        Main method for processing cply files

        Creates Link objects and assigns to member data in Atoms
        """

        import re
        ln = ln.split()
        if (len(ln) >= 4):                
            at = Atom(ln[0], [float(ln[i]) for i in range(1,4)])

            # special 5th column for link tags. parse and store here
            if (len(ln) >= 5):
                tag = ln[4].split("_")
#                    print tag
#                if (tag[0] == "term" or tag[0] == "func" or tag[0] == "rg"):
                if (tag[0] == "term" or tag[0] == "func" or tag[0] == "rg" or tag[0] == "f" or tag[0] == "termfunc" or tag[0] == "termf"):
#                        if (tag[0] == "term"):
#                            at.elmnt = "Al"
                    link = Link(tag[0], cnt, -1)
                    at.link = link
                if (tag[0] == "termcap" or tag[0] == "funccap" or tag[0] == "rgcap" or tag[0] == "fcap"):
                    tag2 = re.split('[()]', tag[3])
                    # tag2[1] now holds which atom this cap is a cap for (1-based) 
#                    print "tag2 ", tag2
                    cap_target = int(tag2[1]) - 1
                    link = Link (tag[0], cnt, cap_target)
                    at.link = link

                # make_ff 
                if( tag[0] == "term" ):
                    at.ctype = "T"
                if( tag[0] == "termcap" ):
                    at.ctype = "X"
                if( tag[0] == "func" ):
                    at.ctype = "F"
                if( tag[0] == "funccap" ):
                    at.ctype = "R"
                    
            self.atoms.append(at)                


    def dist_to_other_structure(self, frag, exclude1 = [], exclude2 = []):
        """
        exclude = list of atoms in frag whose distances are not relevant (b/c they will be deleted
        if we are in midst of a connection process, and these are the "template" hydrogens)
        input in exclude list are indices for Atoms
        """
        dmin = 10000
        for i in range(0, len(self.atoms)):
            if (i not in exclude1):
                for j in range(0, len(frag.atoms)):
                    if (j not in exclude2):
                        d = vec_norm(vec_sub(self.atoms[i].pos, frag.atoms[j].pos))
                        dmin = min(d,dmin)
                        imin = i
                        jmin = j
#        print "dist to other: imin, jmin, dmin = ", imin, jmin, dmin, ",  exclude1, exclude2 = ", exclude1, exclude2
        return dmin


    def atoms_too_close(self):
        """
        like above, but can be faster for quick failure, so I leave it in
        Testing structures that are already connected

        SWS: double check, may not be called
        """
        natoms = len(self.atoms)
        for i in range(0,natoms):
            a1 = self.atoms[i]
            for j in range(i+1,natoms):
                a2 = self.atoms[j]
                d = a1.dist(a2)
                if (d < too_close_dist):
                    print "atoms %d, %d are %f apart" % (i,j,d), "   ", a1.pos, a2.pos 
                    return True
        return False


    def append_fragments(self, frags):
        """
        now list has been checked, just copy it to ourselves
        called after rotation into place and checked new locations are ok
        SWS: double check called 

        Bookingkeeping for indices when addition done

        Args:
            list of fragments (may only be called with one fragment
        """
        import copy
        for i in range(0,len(frags)):
            curlen = len(self.atoms)
            frag = frags[i] 
#            print "Frag: len ", len(frag.atoms)
            for a in frag.atoms:
                newa = copy.copy(a)
                if (newa.link != None):
                    newa.link.shift(curlen)
                self.atoms.append(newa)
#                print a.elmnt, a.pos, a.link

    def update_links(self, idx):
        """
        Keep track of link index book-keeping
        SWS: check Only called on deleting atoms
        """

        for i in range(0,len(self.atoms)):
            if (self.atoms[i].link != None):
                if (idx < self.atoms[i].link.mypos):
                    self.atoms[i].link.mypos -= 1
#                if (idx == self.atoms[i].link.mypos):
#                    print "WARNING: deleted atom that is part of link, expect crash"
                if (idx < self.atoms[i].link.targetpos):
                    self.atoms[i].link.targetpos -= 1
#                if (idx == self.atoms[i].link.targetpos):
#                    print "WARNING: deleted atom that is part of link, expect crash"

    def delete_atom(self, idx):
        """ 
        'Public' called from donoracceptorsystems(?)
        """
#        print "deleting atom index ", idx
        del self.atoms[idx]
        self.update_links(idx)
    
    def delete_link(self, idx):
        """ 
        'Public' called from donoracceptorsystems(?)
        For terminal atoms after connection made
        """
#        print "deleting link for atom index ", idx
        self.atoms[idx].link = None


    def get_normal(self, root_idx = 0):
        """
        Driver for get_normal_from_PCA
        """
        A = [self.atoms[i].pos for i in range(len(self.atoms))]
        v = get_normal_from_PCA(A)
        return v


    def count_funcs(self):
        """
        Count how many functional groups (attachment points eg multiple points)
        can be connected to structure
        by looking at cply and finding ('non' term and cap)
        """
        cnt = 0
        for i in range(0,len(self.atoms)):
            if (self.atoms[i].link != None):
                if self.atoms[i].link.type == "func" or self.atoms[i].link.type == "f":
                    cnt += 1
        return cnt


    def count_funccaps(self):
        """
        Counts cap
        """
        cnt = 0
        for i in range(0,len(self.atoms)):
            if (self.atoms[i].link != None):
                if self.atoms[i].link.type == "funccap" or self.atoms[i].link.type == "fcap":
                    cnt += 1
        return cnt

    def count_funccaps_for(self, target):
        """
        count how many caps there are for particular target func or term
        Note, including both, they can now (9/25/13) be mixed
        """
        cnt = 0
        for i in range(0,len(self.atoms)):
            if (self.atoms[i].link != None):
                if ((self.atoms[i].link.type == "funccap" or self.atoms[i].link.type == "fcap") and self.atoms[i].link.targetpos == target):
                    cnt += 1
                elif (self.atoms[i].link.type == "termcap" and self.atoms[i].link.targetpos == target):
                    cnt += 1
        return cnt



    # make_ff TK
    #   functions start

    def set_q(self,q):
        import numpy
        import sys 
        
        a_cnt = -1 
        for a in self.atoms:
            a_cnt += 1 
            a.q =  q[a_cnt]
            
    def set_resnumb(self,resnumb):
        import numpy
        import sys 
        
        a_cnt = -1 
        for a in self.atoms:
            a_cnt += 1 
            a.resnumb =  resnumb[a_cnt]
            
    def set_restype(self,restype):
        import numpy
        import sys 
        
        a_cnt = -1 
        for a in self.atoms:
            a_cnt += 1 
            a.restype =  restype[a_cnt]
            
    def setall_resnumb(self,resnumb):
        import numpy
        import sys 
        
        for a in self.atoms:
            a.resnumb =  resnumb
            
    def setall_restype(self,restype):
        
        for a in self.atoms:
            a.restype =  restype.strip()
            
    def return_asymb(self):
        import numpy
        
        ASYMB = [] 
        for a in self.atoms:
            ASYMB.append( a.elmnt )
    
        return ASYMB    

    def return_r(self):
        import numpy
        
        R = [] 
        for a in self.atoms:
            r_i =  numpy.array(  [a.pos[0], a.pos[1], a.pos[2]] )
            R.append( r_i )
    
        return R

    def return_ctype(self):
        import numpy
        
        ctype = [] 
        for a in self.atoms:
            ctype.append( a.ctype )
    
        return ctype
    
    
    def return_resnumb(self):
        import numpy
        
        a_cnt = 0
        
        resnumb = [] 
        for a in self.atoms:
            resnumb.append( a.resnumb )
            
        return resnumb    
    
    
    def return_restype(self):
        import numpy
        
        restype = [] 
        for a in self.atoms:
            restype.append( a.restype )
            
        return restype    

    def return_q(self):
        import numpy
        
        q = [] 
        for a in self.atoms:
            q.append( a.q )
            
        return q    

            

class RotationControl(object):
    """
    This class captures options for how to handle the remaining degree of 
    freedom (around the bond direction) left when we slide two fragments into position
    with respect to eachother
    """
    def __init__(self, free_rot_criteria, pattern_str):
        """
        Args:
            free_rot_criteria 
            pattern_str       overrides regular criteria
        """
        self.free_rot_criteria = free_rot_criteria
        self.counter = 0
        self.alternation_pattern = None
        if (pattern_str != None):
            self.alternation_pattern = [float(f) for f in pattern_str.split()]



class Fragment(Structure):
    """
    
    """

    def __init__(self, fname=None, frag=None):
        """
        Constructor:

        Args: 
            fname: filename of a building block cply 
            frag:  existing fragment object (copy constructor)
        """
        self.atoms = []
        if (fname != None):
            self.read(fname)
            self.find_neighbors()
        elif (frag != None):
            Structure.__init__(self, frag)
        else:
            Structure.__init__(self)


    def dump(self):
        """
        For debug
        """
        print "dumping fragment derived from ", self.fname 
        for i in range(0,len(self.atoms)):            
            s = "%d      %s %s  " % (i, str(self.atoms[i].pos), self.atoms[i].elmnt)
            if (self.atoms[i].link != None):
                s += self.atoms[i].link.dump(return_str=True)
            print s


    def get_connectors(self, cap_idx, term_tag, cap_tag):
        """
        this function is written in terms of caps.
        First it finds the cap_idx'th cap of type cap_tag.
        This determines a target atom index that this cap is meant to attach to.
        It returns the term, func, or termfunc index and the cap index corresponding to cap_idx'th
        return index in structure of term_idx'th atom with tag 'term',
        and index of its corresponding 'cap'

        Args:
            cap_idx  : 
            term_tag :
            cap_tag  :
        """

        # import re
        #  print "capidx, term_tag, cap_tag ", cap_idx, term_tag, cap_tag
        ## first find cap_idx'th cap of the right type
        cap_atom_idx = None
        capcnt = 0
        for i in range(0,len(self.atoms)):
#            if (self.atoms[i].link != None):
#                self.atoms[i].link.dump()
            if (self.atoms[i].link != None and self.atoms[i].link.type in cap_tag):
                if (capcnt == cap_idx):
                    # found correct cap
                    cap_atom_idx = i
                    break
                capcnt+=1

        if (cap_atom_idx == None):
            print "cap index %d not found in fragment, expect crash!" % (cap_idx)

        term_atom_idx = self.atoms[cap_atom_idx].link.targetpos

        if (self.atoms[term_atom_idx].link == None):
            raise ValueError, "target term/func/termfunc index does not have a link!"

        if (self.atoms[term_atom_idx].link.type not in term_tag):
            raise ValueError, "target term/fun/termfunc index'th link type (%s) mismatch (%s)" %(str(self.atoms[term_atom_idx].link.type), str(term_tag)) 

        return term_atom_idx, cap_atom_idx




    def get_term(self, term_idx):
        """
        Called from donor....py
        Following connections in cply file
        """
        return self.get_connectors(term_idx, ["term","termfunc","termf"], ["termcap"])

    def get_func(self, term_idx):
        """
        Called from donor....py
        Following connections in cply file
        """
        return self.get_connectors(term_idx, ["func","termfunc"], ["funccap"])

    def get_f(self, term_idx):
        """
        Called from donor....py
        Following connections in cply file
        """
        return self.get_connectors(term_idx, ["f","termf"], ["fcap"])

    def get_rgroup(self, term_idx):
        """
        Called from donor....py
        Following connections in cply file
        """
        return self.get_connectors(term_idx, ["rg"], ["rgcap"])



    def find_neighbors(self, exclude=[]):
        """
        find all neighbors of all atoms
        com ---> "center of mass", will be used to place oxygens at end
        exclude is list of indices that don't count as neighbors (allow substructure search)
        Simple distance search
        """

        self.com = [0,0,0]
        tot_at = 1
        natoms = len(self.atoms)
        for i in range(0,natoms):
            a1 = self.atoms[i]
            if (i not in exclude):
                self.com = vec_add(self.com, a1.pos)
                tot_at += 1
            for j in range(i+1,natoms):
                a2 = self.atoms[j]
                d = a1.dist(a2)
                # print i,j, d
                if (d < bond_dist1 and i not in exclude and j not in exclude):                    
                    a1.nb1.append(j)
                    a2.nb1.append(i)
                    a1.nb.append(j)
                    a2.nb.append(i)
                elif (d < bond_dist2 and i not in exclude and j not in exclude):
                    a1.nb2.append(j)
                    a2.nb2.append(i)                    
                    a1.nb.append(j)
                    a2.nb.append(i)
                # special: also keep track of all neighbors, including spinach
                if (d < bond_dist2):
                    a1.nball.append(j)
                    a2.nball.append(i)

        self.com = [self.com[i]/float(tot_at) for i in range(0,3)]



    def add_coplanar(self, idx, atm):
        """
        given triplet of atoms/indices, find all other atoms on the same plane
        SWS: check not called
        """
        plane_atm = atm
        plane_idx = idx
        # for i in range(0,len(self.atoms)):
        for i in range(0,60):
            if (not i in idx and i < 60):
                if (coplanar(atm, self.atoms[i].pos)):
                    plane_atm.append(self.atoms[i])
                    plane_idx.append(i)
        return plane_idx, plane_atm


    def connected(self, test_idx):
        """
        not full check. enough(?) to just see that
        every atom has at least one neighbor that is in the set as well

        SWS: double check not called
        """
        for i1 in test_idx:
            ok = False
            for n in self.atoms[i1].nb:
                if (n in test_idx):
                    ok = True
            if (not ok):
                return False
        return True


    def bonded(self, i1, i2):
        """
        are two atoms neighbors?
        """
        return (i1 in self.atoms[i2].nb)


    def choose_2pt_free_rotation_angle(self, TR, frag, free_axis, zero_pos, criteria, bond_idx):
        #        criteria = "match_normal"
        #        criteria = "stick_up"
        #        criteria = "max_dist_to_com"
        """
        Args:
            TR        : translation/rotation object
            frag      : fragment object attempting to attach to 'this' fragment object
            free_axis : axis rotation 'frag' around
            zero_pos  : position
            criteria 
                 criteria = "match_normal"
                 criteria = "stick_up"
                 criteria = "max_dist_to_com"
            bond_idx  : small matrix (2x2)
        """
        apply_jitter = True
        # jitter_eps = 30.0 * pi / 180.0  # i.e. 2 degrees out of perfect planarity, so Gaussian notices the degree of freedom
        jitter_eps = 30.0 * math.pi / 180.0  # i.e. 2 degrees out of perfect planarity, so Gaussian notices the degree of freedom
        target_normal = self.get_normal(root_idx=bond_idx[0][0])
        # print "choosing attachment angle around axis ", free_axis, ", criteria = ", criteria, ",  bond_idx =", bond_idx
        # print "my normal = ", target_normal, ", frag normal = ", frag.get_normal(root_idx=bond_idx[0][1]) 
        nj = 40
        maxmindistcom = 0
        maxmindistall = 0
        mindot = 10
        maxdot = 0
        mintheta = None
        testfrag = Fragment(frag=frag)
        natoms = len(testfrag.atoms)

        if (natoms < 3):
            return 0   # nothing to rotate!

        # pick by scan of different angles, pick one by _maximizimng_ minimal distance of atoms in 
        # e.g. rotated indene from center of mass of c60
        for j in range(0,nj):    # discretized angle for 'final' dof  (part of Fragments)
            mindistcom = 1000
            mindistall = 1000
            # theta = 2.0 * pi * float(j) / float(nj)
            theta = 2.0 * math.pi * float(j) / float(nj)
            R2 = get_axis_angle_rotation(free_axis, theta)

            for i in range(0,len(frag.atoms)):       # moving all atoms with trial theta (part of Fragments)
                a = frag.atoms[i]
                ta = apply_transforms(TR, R2, a, zero_pos)
                testfrag.atoms[i].pos = ta
                d = vec_norm(vec_sub(ta, self.com))
                mindistcom = min(d, mindistcom)

            mindistall = self.dist_to_other_structure(testfrag, exclude1 = [bond_idx[1][0]], exclude2 = [bond_idx[0][1]]) #, bond_idx[1][1]])
            test_normal = testfrag.get_normal(root_idx=bond_idx[0][1])
            newdot = (1-abs(vec_dot(target_normal, test_normal)))
            # print "theta = ", theta, ", my normal = ", target_normal, ", test normal = ", test_normal, " newdot = ", newdot 

            if (mintheta == None 
                or (criteria == "max_dist_to_com" and mindistcom > maxmindistcom)                
                or (criteria == "match_normal" and newdot < mindot)
                or (criteria == "stick_up" and newdot > maxdot)
                or (criteria == "max_dist_to_frag" and mindistall > maxmindistall)):
                mintheta = theta
                maxmindistcom = mindistcom
                maxmindistall = mindistall
                mindot = newdot 
                maxdot = newdot
                # print "crit, mindistall", criteria, mindistall
                # testfrag.dump()

        # print "maxmindistall = ", maxmindistall
        # print "mintheta = ", mintheta
        theta = mintheta
        if (apply_jitter and criteria=="match_normal"):
            theta += jitter_eps
        return (theta)



    def place_fragment(self, frag, bond_idx, rot_ctrl):
        """
        Driver for 'several above methods' SWS: Check

        still assuming len(bond_idx) == 2.

        This function handles fairly generically the connection of 2 atoms on one fragment to 2 atoms on another.
        this leaves a degree of freedom.
        For indene case we handle by scanning.  in this case the criteria is keeping
        the placed fragment away from the center of mass.
        generically we will just set the rotation zero.  soon we call an overridable member function for this
        and other such task
        job is transform self. "connector atoms" to those of incoming frag, and subsequentaly all other frag atoms

        e.g.:  bond_idx = [[0,3],[4,7]] means self.atoms[0] connects to frag.atoms[3], and self.4 to frag.7 
        """

        selfpos1 = self.atoms[bond_idx[0][0]].pos
        selfpos2 = self.atoms[bond_idx[1][0]].pos
        fragpos1 = frag.atoms[bond_idx[0][1]].pos
        fragpos2 = frag.atoms[bond_idx[1][1]].pos
        # get transform for given connector points
        TR = get_2pt_transform(selfpos1, selfpos2, fragpos1, fragpos2) 
        # but there is still a degree of freedom:
        free_axis = vec_normalize(vec_sub(selfpos1, selfpos2))

        if (rot_ctrl.free_rot_criteria == "alternate"):
            if (rot_ctrl.counter % 2 == 0):
                free_axis_rot = 0
            else:
                free_axis_rot = pi
        elif (rot_ctrl.alternation_pattern != None):
            # print "alt pat:", rot_ctrl.alternation_pattern
            lpat = len(rot_ctrl.alternation_pattern)
            rot_idx = rot_ctrl.counter % lpat
            free_axis_rot = rot_ctrl.alternation_pattern[rot_idx]*(pi/180.0)
        else:
            free_axis_rot = self.choose_2pt_free_rotation_angle(TR, frag, free_axis, selfpos1, rot_ctrl.free_rot_criteria, bond_idx)
        
        R2 = get_axis_angle_rotation(free_axis, free_axis_rot)
        # now, with optimal orientation of thing we're attaching, copy and transform 
        # to a new fragment
        newfrag = Fragment(frag=frag)
        for i in range(0,len(frag.atoms)):
            a = frag.atoms[i]
            ta = apply_transforms(TR, R2, a, selfpos1)
            e = a.elmnt
            newfrag.atoms[i].pos = ta
        
#        newfrag.dump()
        return newfrag



class FragmentWithHeader(Fragment):
    """
    Header info in cply file (eg  D1(R1,R2,f1,f1)
    """

    def __init__(self, fname=None, frag=None):
        self.atoms = []
        if (fname != None):
            self.read(fname)
            self.find_neighbors()
        elif (frag != None):
            Fragment.__init__(self, frag)
        else:
            Fragment.__init__(self)


    def read(self, fname):
        import re
        self.fname = fname
        # read from piece of Gaussian input file (lines are like xyz: "atom x y z")
        cnt = 0
        first = True
        for ln in fileinput.input(fname):
            if (first):
                # first line is, e.g. "D(R1,R1)"
                # resulting rgroup_tokens is ['R1','R1']
                # we use these in donoracceptorsystems.make_frags() 
                first = False
                ln = ln.strip()
                self.rgroup_spec_str = ln
                p1 = re.split("\(([^)]*)\)", ln)
                p1 = p1[1].split(",")
                self.rgroup_tokens = []
                for p in p1:
                    if (p != ''):
                        self.rgroup_tokens.append(p)
            else:
                self.process_line(ln,cnt)
                cnt += 1
#        self.dump()
#        print self.rgroup_tokens
        self.write_xyz("%s.xyz" % self.fname)
        self.parse_header()


    def parse_header(self):
        """
        gather info from ordered signature list
        eg hdr = ["R1", "R1", "f1", "f2", "f3", "f3", "f1", "f2"]
        unique = ["R1", "f1", "f2", "f3"]
        template = [0 0 1 2 3 3 1 2] = index into unique of each entry in hdr.
        """

        hdr = self.rgroup_tokens
        unique = []
        template = []
        order = {}
        j = 0
        for i in range(len(hdr)):
            tok = hdr[i]
            if (tok not in unique):
                unique.append(tok)
                order[tok] = j
                j+=1
        for i in range(len(hdr)):
            tok = hdr[i]
            template.append(order[tok])
        self.unique = unique
        self.order = order
        self.template = template


if __name__=="__main__":
    # test the "normal from PCA" code.
    atoms = [[1,1,0],[2,1,0],[1,3.4,0],[0,2,0.0]]
    v = get_normal_from_PCA(atoms)
    print "normal for ", atoms
    print "is ", v
