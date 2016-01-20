"""
Class data structures for atomic building blocks 
"""


import sys
import copy
import numpy as np
import math 

import structureContainer as strucC #import StructureContainer
from particles     import Particle #, ParticleContainer
from bonds         import Bond,     BondContainer
from angles        import Angle,    AngleContainer
from dihedrals     import Dihedral, DihedralContainer
from impropers     import Improper, ImproperContainer
from periodictable import periodictable

class Buildingblock(strucC.StructureContainer):

    def __init__(self,name="block", verbose=False):
        """
        Constructor for derived class. The base class constructor is called
        explicitly

        """

        # Base class constructor is called
        strucC.StructureContainer.__init__(self,verbose=False)

        self.verbose = verbose
        self.name = name
        self.composition = []
        self.type = "mol"
        self.connect_id ="term_"
        self.cap_id = "termcap_"
        self.segments = []
        self.func_cnt = 0
        
    def __del__(self):
        """
        Destructor, clears object memory
        """
        #del self

    def find_composition(self):
        """
        Find count of each element
        """
        pt = periodictable()
        self.composition = np.zeros(pt.maxgid)
        
        for pid_i, ptclObj_i  in self.ptclC:
            atomic_number  = ptclObj_i.tagsDict["number"]
            self.composition[atomic_number] += 1 
        
    def setStructureContainer(self, strucC):
        """
        Setter for the structure container.
        """
        self.ptclC  = strucC.ptclC
        self.bondC  = strucC.bondC
        self.angleC  = strucC.angleC
        self.dihC  = strucC.dihC
        self.impC  = strucC.impC         

        for pid_i, ptclObj_i  in self.ptclC:

            # atomic_tag_list   = ["symbol","number","mass","chain","residue","resname","qgroup","fftype","cov_radii","vdw_radii","lmptype","gtype"]
            atomic_tag_list_str = ["symbol","resname","fftype","label"]
            atomic_tag_list_float = ["cov_radii"] #,"vdw_radii"]
            atomic_tag_list_int   = ["number","chain","residue","qgroup","lmptype","ring"]
            for tagid in atomic_tag_list_str :
                if( not isinstance( ptclObj_i.tagsDict[tagid] , str ) ):
                    error_line = "{} tagsDict not set ".format(tagid)
                    sys.exit(error_line)
            for tagid in atomic_tag_list_float :
                if( not isinstance( ptclObj_i.tagsDict[tagid] , float ) ):
                    error_line = "{} tagsDict not set ".format(tagid)
                    sys.exit(error_line)
            for tagid in atomic_tag_list_int :
                if( not isinstance( ptclObj_i.tagsDict[tagid] , int ) ):
                    error_line = "{} tagsDict not set ".format(tagid)
                    sys.exit(error_line)
                    
            # Tags just for building blocks 
            add_dict = ptclObj_i.tagsDict
            add_dict["cplytag"] = ""
            add_dict["linkid"] = ""
            add_dict["link"] = ""
            ptclObj_i.setTagsDict(add_dict)


    def set_cply_tags(self,verbose=False ):
        """
        Set tags for new cply file 
           Use ctype tag and bonding enviroment
        """

        debug = False 

        if( verbose): print "    setting cply tags "
            
        for pid_i, ptclObj_i  in self.ptclC:

            if( debug ):
                print "CHECKING CTYPE ",pid_i,ptclObj_i.tagsDict["linkid"],ptclObj_i.tagsDict["cplytag"]

            # Set terminal attached carbons 
            if( ptclObj_i.tagsDict["linkid"] == "T" and ptclObj_i.tagsDict["number"] == 6 ):
                #
                term_con_cnt = 0 # Terminal connections count
                #
                # Find terminal hydrogen
                #		    
                N_o = cov_nbindx[pid_i]
                N_f = cov_nbindx[pid_i+1] - 1
                for j_indx in range( N_o,N_f+1):
                    pid_j = cov_nblist[j_indx]
                    ptclObj_j = self.ptclC[pid_j]
                    if( ptclObj_j.tagsDict["linkid"] == "X" and ptclObj_j.tagsDict["number"] == 1 ):
                        # Check to see if it has been bonded to another unit 
                        if( ptclObj_j.tagsDict["residue"] ==  ptclObj_j.tagsDict["residue"] ):
                            term_con_cnt += 1
                            ptclObj_j.tagsDict["cplytag"] = "termcap_H(" + str(pid_j) + ")_on_C("+str(pid_i)+")"

                if( term_con_cnt == 1 ):
                    ptclObj_i.tagsDict["cplytag"] = "term_C(" + str(pid_i) + ")"
                if( term_con_cnt > 1 ):
                    print " Number of terminal atoms attached to atom ",pid_i," greater than 1 "
                    sys.exit(" Error in terminal connections ")

            # Set functional attached carbons 
            if( ptclObj_i.tagsDict["linkid"] == "F" and ptclObj_i.tagsDict["number"] == 6 ):
                #
                term_con_cnt = 0 # Terminal connections count
                # 
                # Find function hydrogen
                #
                N_o = cov_nbindx[pid_i]
                N_f = cov_nbindx[pid_i+1] - 1
                for j_indx in range( N_o,N_f+1):
                    pid_j = cov_nblist[j_indx]
                    ptclObj_j = self.ptclC[pid_j]
                    if( ptclObj_j.tagsDict["linkid"] == "R" and ptclObj_j.tagsDict["number"] == 1 ):
                        if( ptclObj_j.tagsDict["residue"] ==  ptclObj_j.tagsDict["residue"] ):
                            term_con_cnt += 1
                            ptclObj_j.tagsDict["cplytag"] = "funccap_H(" + str(pid_j) + ")_on_C("+str(pid_i)+")"

                if( term_con_cnt == 1 ):
                    ptclObj_i.tagsDict["cplytag"] = "func_C(" + str(pid_i) + ")"

                if( term_con_cnt > 1 ):
                    print " Number of functional atoms attached to atom ",pid_i," not equal to 1 "
                    sys.exit(" Error in functional connections ")

        return 	

    def zero_unitq(self,zero_term,zero_func,verbose=False):    
        """
        Sum exsessive charges into carbon atoms 
        """

        for pid_i, ptclObj_i  in self.ptclC:                
            # Find each terminal hydrogen
            if(  ptclObj_i.tagsDict["linkid"]  == "X"  and zero_term  ):
                term = pid_i 
                # Check to be sure hydrogen
                if( ptclObj_i.tagsDict["number"]  != 1 ):
                    print " Non hydrogen used as terminal group atom ",pid_i  
                    sys.exit(" Code unable to process multi atom (nonhyrdogen) terminal group ")
                if( verbose ):
                    print " Terminal atom found ",pid_i," ",ptclObj_i.tagsDict["fftype"]
                #
                # Loop over nieghbors to find attached atom
                #
                N_o = cov_nbindx[pid_i]
                N_f = cov_nbindx[pid_i+1] - 1
                for j_indx in range( N_o,N_f+1):                
                    pid_j = cov_nblist[j_indx]
                    ptclObj_j = self.ptclC[pid_j]

                    term_con_cnt = 0 # Terminal connections count
                    if( ptclObj_j.tagsDict["linkid"].strip() == "T" ):
                        term_con_cnt += 1

                        print "checking ",pid_j,ptclObj_j.tagsDict["linkid"].strip()

                        if( verbose ):
                            print " Terminal connection found ",pid_j," ",ptclObj_j.tagsDict["fftype"]
                        term_con = pid_j

                # Check to be sure multiple atoms not found
                if( term_con_cnt < 1 ):
                    print " No terminal connections found "
                    sys.exit(" Error in terminal connections ")

                if( term_con_cnt > 1 ):
                    print " Multiple terminal connections found "
                    sys.exit(" Error in terminal connections ")

                # Sum charges into base monomer unit
                self.ptclC[term_con].charge = self.ptclC[term_con].charge  + self.ptclC[term].charge 
                self.ptclC[term].charge   = 0.0

            # Find each functional hydrogen
            if( ptclObj_i.tagsDict["linkid"] == "R"  and zero_func ):
                term = pid_i 
                # Check to be sure hydrogen
                if(  ptclObj_i.tagsDict["number"] != 1 ):
                    print " Non hydrogen used as functional group "
                    sys.exit(" Code unable to process multi atom (nonhyrdogen) functional group ")
                if( verbose ):
                    print " Functional atom found ",pid_i," ",ptclObj_i.tagsDict["fftype"]
                #
                # Loop over nieghbors to find attached atom
                #
                N_o = cov_nbindx[pid_i]
                N_f = cov_nbindx[pid_i+1] - 1
                for j_indx in range( N_o,N_f+1):
                    pid_j = cov_nblist[j_indx]
                    ptclObj_j = self.ptclC[pid_j]
                    term_con_cnt = 0 # Terminal connections count
                    if( ptclObj_j.tagsDict["linkid"].strip() == "F" ):
                        term_con_cnt += 1
                        if( verbose ):
                            print " functional connection found ",pid_j," ",ptclObj_j.tagsDict["fftype"]
                        term_con = pid_j
                # Check to be sure multiple atoms not found
                if( term_con_cnt  < 1 ):
                    print " No functional connections found "
                    sys.exit(" Error in functional connections ")
                # Check to be sure multiple atoms not found
                if( term_con_cnt > 1 ):
                    print " Multiple functional connections found "
                    sys.exit(" Error in functional connections ")
                # Sum charges into base monomer unit
                self.ptclC[term_con].charge = self.ptclC[term_con].charge  + self.ptclC[term].charge 
                self.ptclC[term].charge   = 0.0

        return



    def write_cply(self,cply_file, write_ff=False,write_bonds=False,write_latvec=True,write_segments=False):
        """
        Write cply file
        """

        F = open(cply_file,'w')
        cply_line = "# atomic_symb ,float(r_i[0]), float(r_i[1]),float(r_i[2]),label,fftype,ptclObj.mass,charge,qgroup,ring,residue,resname, chain, cplytag"
        F.write(cply_line)
        for pid, ptclObj  in self.ptclC:
            r_i = ptclObj.position
            atomic_symb = ptclObj.tagsDict['symbol']
            cplytag = ptclObj.tagsDict["cplytag"]
            if(write_ff ):
                charge = ptclObj.charge
                residue = ptclObj.tagsDict["residue"]
                resname = ptclObj.tagsDict["resname"]
                label = ptclObj.tagsDict["label"]
                fftype = ptclObj.tagsDict["fftype"]
                qgroup = ptclObj.tagsDict["qgroup"]
                chain = ptclObj.tagsDict["chain"]
                ring = ptclObj.tagsDict["ring"]
                cply_line =  "\n %5s %16.8f %16.8f %16.8f %s %s %12.8f %12.8f  %d %d %d %s  %d %s "%(atomic_symb ,float(r_i[0]), float(r_i[1]),float(r_i[2]),label,fftype,ptclObj.mass,charge,qgroup,ring,residue,resname, chain, cplytag )
            else:
                cply_line =  "\n %5s %16.8f %16.8f %16.8f %s "%(atomic_symb ,float(r_i[0]), float(r_i[1]),float(r_i[2]), cplytag )
            F.write(cply_line)

        if( write_latvec ):
            F.write("\n lattice")
            for d in range(3):
                cply_line = "\n %16.8f %16.8f %16.8f "%( self.latvec[d][0],self.latvec[d][1],self.latvec[d][2])
                F.write(cply_line)
            

        if( write_bonds ):
            for b_o, bondObj_o  in self.bondC:
                pid_i = bondObj_o.pgid1
                pid_j = bondObj_o.pgid2
                F.write("\n  bond %d %d  "%(pid_i,pid_j))

        if( write_segments ):

            for segment_i in self.segments :
                F.write("\n segment {} ".format(segment_i["tag"]))
                F.write("\n unit {} ".format( segment_i["unit"] ))
                func_line = "\n func "
                for func_i in segment_i["func"]:
                    func_line += " {} ".format( func_i )
                F.write(func_line)
                F.write("\n end")
            
        F.close()

        return 

    def read_cply(self,cply_file,debug = False):
        """
        Read cply file
        """
        
        # Load periodic table 
        pt = periodictable()

        read_lattice = False
        read_segment = False 

        with open(cply_file) as f:
            for line in f:
                col = line.split()
                if( read_lattice ):
                    self.latvec[lv_cnt][0] = float( col[0] )
                    self.latvec[lv_cnt][1] = float( col[1] )
                    self.latvec[lv_cnt][2] = float( col[2] )
                    if( lv_cnt == 2 ):
                        read_lattice = False 
                    lv_cnt += 1
                elif( read_segment ):
                    if( debug ):
                        print " Readin degment line ",col
                    
                    if( str(col[0]) == 'unit' ):
                        if( len(col) > 2 ):
                            error_line = "Each segment can only have 1 unit "
                            error_line += "\n {} in unit line will be ignored ".format(str(col[2::]))
                            print error_line
                        segment_i['unit'] =  col[1] 
                    if( str(col[0]) == 'func' ):
                        segment_i['func'] =  col[1::] 
                    if( str(col[0]) == 'end' ):
                        read_segment = False 
                                        
                elif( len(col) >= 4  and col[0] != "bond" and col[0] != "#" ):

                    pt_i = Particle()
                    pt_i.type = str(col[0])
                    pt_i.position = [ float(col[1]),float(col[2]),float(col[3]) ]
                    
                    add_dict = pt_i.tagsDict
                    el = pt.getelementWithSymbol(str(col[0]))
                    add_dict["symbol"] = str(col[0])
                    add_dict["number"] = el.number
                    add_dict["mass"] = el.mass

                    add_dict["cplytag"] = ""
                    add_dict["linkid"] = ""
                    add_dict["link"] = ""
                    pt_i.mass = el.mass
                    add_dict["chain"] = 1
                    add_dict["ring"] = 0
                    add_dict["residue"] = 1
                    add_dict["resname"] = "RES"
                    add_dict["qgroup"] = 1
                    add_dict["fftype"] = "??"
                    add_dict["label"] =  el.symbol  
                    add_dict["cov_radii"] = el.cov_radii
                    add_dict["vdw_radii"] = el.vdw_radii
                    add_dict["lmptype"] = -1
                    if (len(col) >= 14 ):
                        add_dict["label"] = str(col[4])
                        add_dict["fftype"] = str(col[5])
                        add_dict["mass"] = float(col[6])                        
                        pt_i.mass = float(col[6])                   
                        pt_i.charge = float(col[7])
                        add_dict["qgroup"] = int(col[8])                        
                        add_dict["ring"] = int(col[9])                        
                        add_dict["residue"] = int(col[10])
                        add_dict["resname"] = str(col[11])                        
                        add_dict["chain"] = int(col[12])                        
                        add_dict["cplytag"] = str(col[13])
                    elif (len(col) == 13 ):
                        add_dict["label"] = str(col[4])
                        add_dict["fftype"] = str(col[5])
                        add_dict["mass"] = float(col[6])                        
                        pt_i.mass = float(col[6])                   
                        pt_i.charge = float(col[7])
                        add_dict["qgroup"] = int(col[8])                        
                        add_dict["ring"] = int(col[9])                        
                        add_dict["residue"] = int(col[10])
                        add_dict["resname"] = str(col[11])                        
                        add_dict["chain"] = int(col[12])                        
                    elif (len(col) == 8 ):
                        add_dict["residue"] = int(col[5])
                        add_dict["resname"] = str(col[6])                        
                        add_dict["cplytag"] = str(col[7])
                    elif (len(col) == 7 ):
                        pt_i.charge = float(col[4])
                        add_dict["residue"] = int(col[5])
                        add_dict["resname"] = str(col[6])                                                
                    elif (len(col) == 5 ):
                        add_dict["cplytag"] = str(col[4])

                    pt_i.setTagsDict(add_dict)                    
                    self.ptclC.put(pt_i)

                    # # # print "debug pt_i ",len(col),pt_i
                    
                elif(len(col) >= 3 ):
                    if( col[0] == "bond"):
                        b_i = int(col[1])
                        b_j = int(col[2])
                        bnd = Bond(b_i,b_j)
                        #print "process_line bond line ",col
                        self.bondC.put(bnd)
                # Key word search 
                if( len(col) > 0 ):
                    if( str(col[0]) == 'lattice' ):
                        read_lattice = True
                        lv_cnt = 0
                    if( str(col[0]) == 'segment' ):
                        read_segment = True 
                        segment_i = {} #segment(str(col[1]))
                                       # segment_i = segment()
                        segment_i['tag'] = str(col[1])
                        segment_i['segment'] = {}
                        self.segments.append(segment_i)
        if( debug ):
            print segments
            sys.exit("debug in read_cply is True ")
        
    def set_type(self):
        """
        Count the number of connection points
        and set the type to unit or term depending 
        """
        self.cnt_connectors()
        if( self.connect_cnt == 1 ):
            self.type = "term"
        elif( self.connect_cnt == 2 ):
            self.type = "unit"
        else:
            error_line = " Unknow type with {} connectors ".format(self.connect_cnt)
            print error_line
                        
    def cnt_connectors(self):
        """
        Count the number of connection points
        and set the type to unit or term depending 
        """
        self.connect_cnt = 0 
        self.func_connect_cnt = 0 
        for pid_i, ptclObj_i  in self.ptclC:    
            if( ptclObj_i.tagsDict["cplytag"][:5] == "term_" ):
                self.connect_cnt += 1
            if( ptclObj_i.tagsDict["cplytag"][:5] == "func_" ):
                self.func_connect_cnt += 1

    def get_connectors(self,connect_indx, verbose=False ):
        """
        Get connector index numbers for 
        """
        if( verbose ):
            print " Checking for terminal {}".format(connect_indx)

        connect_cnt = 0
        self.pid_term = -1 
        self.pid_termcap = -1 
        for pid_i, ptclObj_i  in self.ptclC:    
            if( ptclObj_i.tagsDict["cplytag"][:5] == self.connect_id ):
                if( verbose ):
                    print connect_cnt,"    {}  ",self.connect_id, ptclObj_i.tagsDict['symbol'],ptclObj_i.tagsDict["cplytag"],ptclObj_i.tagsDict["linkid"],ptclObj_i.tagsDict["link"]
                
                if( connect_cnt == connect_indx ):
                    self.pid_term = pid_i
                    N_o = self.bonded_nbindx[pid_i ]
                    N_f = self.bonded_nbindx[pid_i+ 1 ] - 1
                    for indx_j in range( N_o,N_f+1):
                        pid_j = self.bonded_nblist[indx_j]
                        if( self.ptclC[pid_j].tagsDict["cplytag"][:8] == self.cap_id ):
                            if( verbose ):
                                print "    termcap ", ptclObj_i.tagsDict['symbol'],ptclObj_i.tagsDict["cplytag"],ptclObj_i.tagsDict["linkid"],ptclObj_i.tagsDict["link"]
                            self.pid_termcap = pid_j

                    
                connect_cnt += 1
                                
        if( self.pid_term < 0 ):
            for pid_i, ptclObj_i  in self.ptclC:
                print ptclObj_i.tagsDict['symbol'],ptclObj_i.tagsDict["cplytag"]
            error_line = " No terminal particle found for connection {} on atom {} ".format(connect_indx,self.pid_term)
            append = False 
            self.ptclC.write_xmol("error.xyz",error_line,append)
            sys.exit(error_line)
        if( self.pid_termcap < 0 ):
            error_line = " No terminal cap particle found for connection {} ".format(connect_indx)
            append = False 
            self.ptclC.write_xmol("error.xyz",error_line,append)
            pid_i = self.pid_term
            N_o = self.bonded_nbindx[pid_i ]
            N_f = self.bonded_nbindx[pid_i+ 1 ] - 1
            print "i = ",pid_i
            for indx_j in range( N_o,N_f+1):
                pid_j = self.bonded_nblist[indx_j]
                print " nieghbor of i ",pid_j,self.ptclC[pid_j].tagsDict["cplytag"][:8]
            
            sys.exit(error_line)

        
        
        return 
        

    def align_termcaps(self ,verbose=False, debug = False ):
        """
        Align term caps along the x-axis
        """

        def calc_cos(v_i,v_j):
            """
            Calculate the cos theta between two vectors 
            """
            return np.dot(v_i/np.linalg.norm(v_i),v_j/np.linalg.norm(v_j))

        self.cnt_connectors()
        if( self.connect_cnt >= 2 ):
            # Get position of first term cap particle 
            self.get_connectors(0 )
            r_i = np.array(self.ptclC[self.pid_termcap].position )
            # shift origin to r_i
            self.ptclC.shift(-1.0*r_i)
            # Get position of second term cap particle
            self.get_connectors( self.connect_cnt -1 )
            r_j = np.array(self.ptclC[self.pid_termcap].position )
            r_xy = np.array([r_j[0],r_j[1],0.0])
            unit_x = np.array( [1.0,0.0,0.0])
            if( debug ):
                print "r_j",r_j
                print "",r_xy[1]
                print calc_cos(unit_x,r_xy)
                print np.rad2deg(np.arccos(calc_cos(unit_x,r_xy)))
                print "r_xy",r_xy

            if( r_xy[1] > 0.0 ):
                # if in the 1st or 2nd quadrents rotate  clockwise 
                self.rotate_xy(np.arccos(calc_cos(unit_x,r_xy)),direction="clockwise",verbose=False)
            elif( r_xy[1] <= 0.0 ):
                # if in the 3rd or 4th quadrents rotate counter clockwise 
                # of if 180 rotation 
                self.rotate_xy(np.arccos(calc_cos(unit_x,r_xy)),direction="counterclockwise",verbose=False)


            if( debug ):
                print "after xy rot self.ptclC[self.pid_termcap].position",self.ptclC[self.pid_termcap].position
            r_j = np.array(self.ptclC[self.pid_termcap].position )
            r_xz = np.array([r_j[0],0.0,r_j[2]])
            if( debug ):
                print "r_xz",r_xz

            if( r_xz[2] > 0.0 ):
                # if in the 1st or 2nd quadrents rotate  clockwise 
                self.rotate_xz(np.arccos(calc_cos(unit_x,r_xz)),direction="clockwise",verbose=False)
            elif( r_xz[2] <= 0.0 ):
                # if in the 3rd or 4th quadrents rotate counter clockwise
                # of if 180 rotation 
                self.rotate_xz(np.arccos(calc_cos(unit_x,r_xz)),direction="counterclockwise",verbose=False)


            if( debug ):
                print "after xz rot self.ptclC[self.pid_termcap].position",self.ptclC[self.pid_termcap].position

            # Rotat term atom so first term-termcap bond is in the xz plane
            self.get_connectors(0 )
            r_i = np.array(self.ptclC[self.pid_term].position )
            r_zy = np.array([0.0,r_i[1],r_i[2]])
            unit_z = np.array( [0.0,0.0,1.0])
            if( debug ):
                print "r_zy",r_zy
            if( r_zy[1] > 0.0 ):
                # if in the 1st or 2nd quadrents rotate  clockwise 
                self.rotate_yz(np.arccos(calc_cos(unit_z,r_zy)),direction="counterclockwise",verbose=False)
            elif( r_zy[1] <= 0.0 ):
                # if in the 3rd or 4th quadrents rotate counter clockwise 
                # of if 180 rotation 
                self.rotate_yz(np.arccos(calc_cos(unit_z,r_zy)),direction="clockwise",verbose=False)
            # if  r_xy[1] == 0.0 no rotation is needed

            if( debug ):
                print "after zy rot self.ptclC[self.pid_term].position",self.ptclC[self.pid_term].position


            # Shift so second term is at origin
            #   as this will be adding point for the next buildingblock unit  
            #connect_indx = 1 
            #self.get_connectors(connect_indx )
            #self.ptclC.shift(-1.0*np.array(self.ptclC[self.pid_term].position))
        
        return 

    def rm_cplytags(self,verbose=False,debug=False):
        """
        Set all  cplytag to empty string 
        """
        for pid, ptclObj  in self.ptclC :
            ptclObj.tagsDict["cplytag"] = ""

    def get_max(self,verbose=False,debug=False):
        """
        get max chain, residue and charge group number 
        """
        self.max_chain = 0 
        self.max_residue = 0 
        self.max_qgroup = 0 
        self.max_ring = 0 
        for pid, ptclObj  in self.ptclC :
            if( ptclObj.tagsDict["chain"] > self.max_chain ): self.max_chain =  ptclObj.tagsDict["chain"] 
            if( ptclObj.tagsDict["residue"] > self.max_residue ): self.max_residue =  ptclObj.tagsDict["residue"] 
            if( ptclObj.tagsDict["qgroup"] > self.max_qgroup ): self.max_qgroup =  ptclObj.tagsDict["qgroup"]
            if( ptclObj.tagsDict["ring"] > self.max_ring ): self.max_ring =  ptclObj.tagsDict["ring"]
                
    def get_min(self,verbose=False,debug=False):
        """
        get max chain, residue and charge group number 
        """
        self.min_chain = 100000000
        self.min_residue = 100000000
        self.min_qgroup = 100000000
        self.min_ring = 100000000
        for pid, ptclObj  in self.ptclC :
            if( ptclObj.tagsDict["chain"] < self.min_chain ): self.min_chain =  ptclObj.tagsDict["chain"] 
            if( ptclObj.tagsDict["residue"] < self.min_residue ): self.min_residue =  ptclObj.tagsDict["residue"] 
            if( ptclObj.tagsDict["qgroup"] < self.min_qgroup ): self.min_qgroup =  ptclObj.tagsDict["qgroup"]
            if( ptclObj.tagsDict["ring"] < self.min_ring ): self.min_ring =  ptclObj.tagsDict["ring"]
                
    def get_maxmin(self,verbose=False,debug=False):
        self.get_max()
        self.get_min()
                                
    def shift_tag(self,tag,tag_min):
        """
        shift tag by an number  
        """
        for pid, ptclObj  in self.ptclC :
            if( ptclObj.tagsDict[tag] > 0 ):
                ptclObj.tagsDict[tag] += tag_min

    def add_chain(self,max_ref_chain ):
        """
        Set minimum chain number to references

        
        """
        for pid, ptclObj in self.ptclC :
            ptclObj.tagsDict["chain"] +=   max_ref_chain  - self.max_chain 

    def set_tag(self,tag,tag_value):
        """
        set tag to value 
        """
        for pid, ptclObj  in self.ptclC :
            ptclObj.tagsDict[tag] = tag_value
                
    def align_termbond(self,connector_i,origin="term",verbose=False,debug=False):
        """
        Align terminal bond along the x-axis

        origin
           term - sets the term particle to the origin with the term cap along the positive x-axis
           termcap - sets the term cap particle to the origin with the term along the positive x-axis

              C         H
           /     \    /
        H          C 
        """

        def calc_cos(v_i,v_j):
            """
            Calculate the cos theta between two vectors 
            """
            return np.dot(v_i/np.linalg.norm(v_i),v_j/np.linalg.norm(v_j))

        # Get position of  term cap particle 
        self.get_connectors(connector_i)
        if(origin == "termcap" ):
            pid_i = self.pid_termcap
            pid_j = self.pid_term
        elif(origin == "term"):
            pid_i = self.pid_term
            pid_j = self.pid_termcap
        else:
            error_line = "!!! Error in Buildingblock.align_termbond() !!! \n"
            error_line += "!!! Unknow origin selected {} please select term or termcap !!!".format(origin)
            sys.exit(error_line)
        
        r_i = np.array(self.ptclC[pid_i].position )            
        # shift origin to r_i
        self.ptclC.shift(-1.0*r_i)
        r_j = np.array(self.ptclC[pid_j].position )
        # Get position of second  particle
        r_xy = np.array([r_j[0],r_j[1],0.0])
        unit_x = np.array( [1.0,0.0,0.0])
        if( debug ):
            print "r)j",r_j
            print "",r_xy[1]
            print np.rad2deg(np.arccos(calc_cos(unit_x,r_xy)))

        
        if( r_xy[1] > 0.0 ):
            # if in the 1st or 2nd quadrents rotate  clockwise 
            self.rotate_xy(np.arccos(calc_cos(unit_x,r_xy)),direction="clockwise",verbose=False)
        elif( r_xy[1] <= 0.0 ):
            # if in the 3rd or 4th quadrents rotate counter clockwise 
            # of if 180 rotation 
            self.rotate_xy(np.arccos(calc_cos(unit_x,r_xy)),direction="counterclockwise",verbose=False)
        # if  r_xy[1] == 0.0 no rotation is needed

        if( debug ):
            print "after xy rot self.ptclC[pid_j].position ",self.ptclC[pid_j].position 
        r_j = np.array(self.ptclC[pid_j].position  )
        r_xz = np.array([r_j[0],0.0,r_j[2]])
        if( debug ):
            print "r_xz",r_xz
        
        if( r_xz[2] > 0.0 ):
            # if in the 1st or 2nd quadrents rotate  clockwise 
            self.rotate_xz(np.arccos(calc_cos(unit_x,r_xz)),direction="clockwise",verbose=False)
        elif( r_xz[2] <= 0.0 ):
            # if in the 3rd or 4th quadrents rotate counter clockwise 
            # of if 180 rotation 
            self.rotate_xz(np.arccos(calc_cos(unit_x,r_xz)),direction="counterclockwise",verbose=False)
        # if  r_xy[1] == 0.0 no rotation is needed


        if( debug ):
            print "after xz rot self.ptclC[pid_j].position ",self.ptclC[pid_j].position 

        return 


    def __iadd__(self, other):
        """
        'Magic' method to implement the '+=' operator (eg struc1 += struc2)
        
        Compare global IDs of particles and reassign globalIDs for particle
        container using the max ID between the two lists. Tracks these changes
        for all other (bond, angle, dihedral) containers that reference particleIDs
        """

        self.verbose = False 
        
        # Empty container checks
        if len(other) == 0:             # If struc2 is empty (has no particles)
            return self                 #   simply return unchanged current container
        if len(self) == 0:              # If struc1 (this struc) is empty (has no particles)
            return copy.deepcopy(other) #    full copy other and return

        #
        #  Special method for connecting building blocks
        #
        # Count number of connectors
        connect_along_bond = False
        if( self.type != "mol" and other.type != "mol"):            
            # IF not complete molecules that will just be added
            #   connect along bonds
            connect_along_bond = True 

            self.cnt_connectors()
            other.cnt_connectors()
            if( self.verbose ):
                log_line =  "\n building block 1 type {}  ".format(self.type )
                log_line += "\n   building block 1 has {} connections ".format(self.connect_cnt )
                log_line += "\n   building block 1 has {} func connections ".format(self.func_connect_cnt )
                log_line += "\n building block 2 type {}  ".format(other.type )
                log_line += "\n   building block 2 has {} connections ".format(other.connect_cnt )
                log_line += "\n   building block 2 has {} func connections ".format(other.func_connect_cnt )
                print log_line

            self.connect_id = "term_"
            self.cap_id = "termcap_"
            other.connect_id = "term_"
            other.cap_id = "termcap_"
            if( self.type == "unit" and other.type == "unit"):            
                self.connectionpoint = 1 
                other.connectionpoint = 0
                new_termpoint = 1
            elif( self.type == "term" and other.type == "unit"):            
                # elif( self.connect_cnt == 1 and other.connect_cnt  == 2 ):
                self.connectionpoint = 0
                other.connectionpoint = 0
                new_termpoint = 0
            elif( self.type == "term" and other.type == "term"):            
                #        elif( self.connect_cnt == 1 and other.connect_cnt  == 1 ):
                self.connectionpoint = 0
                other.connectionpoint = 0
                new_termpoint = -1
            elif( self.type == "unit" and other.type == "term"):            
                #elif( self.connect_cnt == 2 and other.connect_cnt  == 1 ):
                self.connectionpoint = 1
                other.connectionpoint = 0            
                new_termpoint = -1
            elif( self.type == "unit" and other.type == "func"):            
                #elif( self.connect_cnt == 2 and other.connect_cnt  == 1 ):
                self.connect_id = "func_"            
                self.cap_id = "funccap_"
                # this is zero since each functional possition will be removed once a functional group is added 
                self.connectionpoint = 0            
                other.connect_id = "term_"
                other.connectionpoint = 0            
                new_termpoint = -1
            else:
                error_line = " Unknow connection configuration {} - {}  \n".format(self.type, other.type )
                error_line += " building block 1 has {} connnections \n".format(self.connect_cnt)
                error_line += " building block 2 has {} connnections \n".format(other.connect_cnt)
                sys.exit(error_line)
            
            if( self.verbose ):
                log_line = "\n Adding connection along bond "
                log_line += "\n  Connecting block 2 {} at point {} ".format(other.connect_id,other.connectionpoint )
                log_line += "\n  to  block 1 {} at point {} ".format(self.connect_id,self.connectionpoint )
                print log_line

            # Shift second buildingblock so first connector is bonded to the second connector of  first buildingblock
            if( self.verbose ): print " Shift  along bonds "
            self.align_termbond(self.connectionpoint,origin="term")
            other.align_termbond(other.connectionpoint,origin="termcap")

            if( self.verbose ): print " Find terminal atoms  "
            self.get_connectors(self.connectionpoint)
            other.get_connectors(other.connectionpoint)
            bond_length = np.array([self.ptclC[self.pid_term].radii + other.ptclC[other.pid_term].radii,0.0,0.0])

            # Shift second building block to have the correct bond length 
            other.ptclC.shift((bond_length-1.0*np.array(other.ptclC[other.pid_term].position)))

            # Delete termcaps
            self.bondC.deletepid(self.pid_termcap)
            other.bondC.deletepid(other.pid_termcap)
            del self.ptclC[self.pid_termcap]
            del other.ptclC[other.pid_termcap]

            # Reset term tags to null
            self.ptclC[self.pid_term].tagsDict["cplytag"] = ""
            other.ptclC[other.pid_term].tagsDict["cplytag"] = ""

            # If an intern ring connection is made change the fftype
            if( self.ptclC[self.pid_term].tagsDict["ring"] > 0 and other.ptclC[other.pid_term].tagsDict["ring"] > 0 ):
                #print " ring self ", self.ptclC[self.pid_term].tagsDict["ring"]
                #print " ring other ", other.ptclC[other.pid_term].tagsDict["ring"]
                
                self.ptclC[self.pid_term].tagsDict["fftype"] = 'C!'
                other.ptclC[other.pid_term].tagsDict["fftype"] = 'C!'
                
                # sys.exit(" ring test ")
            
        idFromToDict = dict()  # Need to keep track of all ptcl ID changes at once
                               # {fromID1:toID1, fromID2:toID2...}
                               # eg {1:3, 3:5, 2:20...}
        
        bondC  = BondContainer()              # Local bond container copy so ptclIDs
        angleC = AngleContainer()             # Local angle container copy so ptclIDs
        dihC   = DihedralContainer()          # Local dihedral container copy so ptclIDs
        
        bondC  = copy.deepcopy(other.bondC)   #  inside can be changed (for adding below)
        angleC = copy.deepcopy(other.angleC)  #  inside can be changed (for adding below)
        dihC   = copy.deepcopy(other.dihC)    #  inside can be changed (for adding below)
        impC   = copy.deepcopy(other.impC)    #  inside can be changed (for adding below)
        
        keys1 = self.ptclC.particles.keys()    # global IDs of particles in this object
        keys2 = other.ptclC.particles.keys()   # global IDs in object being added
        self.ptclC.maxgid= max(keys1 + keys2)  # find max globalID in keys, set this object maxID

        self.get_max()
        other.get_max()
        # Update chain, residue and charge group number
        if( self.max_chain == other.max_chain ):

            other.shift_tag("qgroup",self.max_qgroup )
            other.shift_tag("residue",self.max_residue )
            other.shift_tag("ring",self.max_ring )
        else:
            other.add_chain(self.max_chain )
        
        for ptclkey2 in other.ptclC.particles:
            self.ptclC.put(other.ptclC.particles[ptclkey2]) # Pushes ptcl to this struc's ptcl container
            fromPtclID = ptclkey2                           # Track IDs from--->to
            toPtclID   = self.ptclC.maxgid                  #  --> toID (to is the maxid of this ptclC)

            idFromToDict[fromPtclID]=toPtclID               # Store ID changes

        if( connect_along_bond ):
            # Add inter building block bond
            if( self.verbose ): print " Add inter building block bond "
            bb_bb = Bond(self.pid_term,idFromToDict[other.pid_term])
            self.bondC.put(bb_bb)

        bondC.replacePtclIDs(idFromToDict)  # Use tracked list of ID changes
        self.bondC += bondC                 # Now add bondC with 'corrected' IDs

        angleC.replacePtclIDs(idFromToDict) # Now add angleC with 'corrected' IDs
        self.angleC += angleC               # Use tracked list of ID changes

        dihC.replacePtclIDs(idFromToDict)   # Use tracked list of ID changes
        self.dihC += dihC                   # Now add dihC with 'corrected' IDs

        impC.replacePtclIDs(idFromToDict)   # Use tracked list of ID changes
        self.impC += impC                   # Now add impC with 'corrected' IDs

        
        if( connect_along_bond ):
            self.compressPtclIDs()
            self.bondC_nblist()
            if ( new_termpoint >= 0 ):

                if( self.verbose ): print " Add new cplytag to  "
                    
                self.get_connectors(new_termpoint)
                self.ptclC[self.pid_term].tagsDict["cplytag"] = "term_C({})".format(self.pid_term)
                self.ptclC[self.pid_termcap].tagsDict["cplytag"] = "termcap_H({})_on_C({})".format(self.pid_termcap,self.pid_term)


        # Append segments 
        for seg_i in other.segments:
            self.segments.append(seg_i)
        
        return self

