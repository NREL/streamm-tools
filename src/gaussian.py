#! /usr/bin/env python

"""
Subroutines for reading and writing gaussian related files
"""

# Dr. Travis Kemper
# NREL
# Initial Date 06/20/2013
# travis.kemper@nrel.gov

import numpy as np
    

from periodictable import periodictable
from particles  import Particle


def read_fchk(strucC,fchk_file):
    """
    Read in structure information from gaussian fchk file

    Args:
        strucC (str) 
        fchk_file (str) gaussian fchk file

    """



    # Check to see if a previous read has occured
    pt_update = False
    if( len(strucC.ptclC) > 0 ):
        pt_update = True

        
    # Energy conversion
    # http://physics.nist.gov/cgi-bin/cuu/Value?threv
    HtoeV = 27.211385

    bohr2angstrom = 0.5291772086

    F = open(fchk_file,'r')
    Lines = F.readlines()
    F.close()


    # Load periodic table 
    pt = periodictable()

    read_r = False
    read_eln = False
    read_esp = False
    for line in Lines :
        col = line.split()

        if( read_r ):

            if (  col[0] == "Force" and col[1] == "Field" ):
                read_r = False
                p_i = 0 
                for p_indx in range(NA):
                    #print atom_i ,atom_i*3,atom_i*3+2,R_all[atom_i*3:atom_i*3+3]
                    vec_r_i =  R_all[p_indx*3:p_indx*3+3]
                    p_i += 1
                    strucC.ptclC[p_i].position = vec_r_i
                        
            else:
                for r_i in  map(float,col) :
                    R_all.append( r_i*bohr2angstrom )

        if( read_eln ):
            if ( eln_p_cnt == NA ):
                read_eln = False
            else:                    
                for eln_i in  map(int,col):
                    eln_p_cnt += 1
                    if( pt_update ):
                        
                        if( eln_i != strucC.ptclC[eln_p_cnt].tagsDict["number"] ):
                            error_line =  "  Particle ",eln_p_cnt, strucC.ptclC[eln_p_cnt].tagsDict["number"]," != ",eln_i
                            sys.exit(error_line)
                    else:
                        el = pt.getelementWithNumber(eln_i)
                        strucC.ptclC[eln_p_cnt].tagsDict["number"] = eln_i
                        
                        strucC.ptclC[eln_p_cnt].tagsDict["gtype"] = el.symbol
                        strucC.ptclC[eln_p_cnt].tagsDict["fftype"] = "??" 
                        strucC.ptclC[eln_p_cnt].tagsDict["resname"] = "??"
                        strucC.ptclC[eln_p_cnt].tagsDict["residue"] = 0
                        strucC.ptclC[eln_p_cnt].tagsDict["qgroup"] = 0
                        strucC.ptclC[eln_p_cnt].tagsDict["chain"] = 0
                        strucC.ptclC[eln_p_cnt].tagsDict["symbol"] = el.symbol
                        strucC.ptclC[eln_p_cnt].tagsDict["number"] = el.number
                        strucC.ptclC[eln_p_cnt].tagsDict["cov_radii"] = el.cov_radii
                        strucC.ptclC[eln_p_cnt].tagsDict["vdw_radi"] = el.vdw_radii
                        
                        
        if( read_esp ):
            if ( esp_p_cnt == NA ):
                read_esp = False
            else:

                for q_i in  map(float,col):
                    esp_p_cnt += 1
                    p_indx = esp_p_cnt -1 
                    strucC.ptclC[esp_p_cnt].charge = q_i


        if( len(col) > 2 ):
            if( col[0] == "Total" and col[1] == "Energy" ):
                TOTAL_ENERGY = float( col[3] )*HtoeV

        if( len(col) == 5 ):
            if( col[0] == "Number" and col[1] == "of"  and col[2] == "atoms" ):
                NA = int(col[4])
                if( pt_update ):
                    if( NA != len(strucC.ptclC) ):
                        print " json file contains %d atoms and fchk file contains %d "%(len(strucC.ptclC),NA)
                        sys.exit("inconsistent files ")
                else:
                    # Create particles to be updated 
                    for pid_i in range(1,NA+1):
                        pt_i = Particle( [0.0,0.0,0.0] )
                        strucC.ptclC.put(pt_i)

        if( len(col) == 6 ):
            if( col[0] == "Current" and col[1] == "cartesian"  and col[2] == "coordinates" ):
                read_r = True
                R_all = []


        if( len(col) > 2  ):
            if( col[0] == "Atomic" and col[1] == "numbers"   ):
                read_eln = True
                eln_p_cnt = 0

        if( len(col) > 2  ):
            if( col[0] == "ESP" and col[1] == "Charges"   ):
                read_esp = True
                esp_p_cnt = 0

    return strucC,TOTAL_ENERGY

    
def get_dih_id( zmatrix ):
    """
    Find dihedral id's atoms and values from zmatrix
    """
    
    import sys
    
    DIH_ID = []
    DIH_VAL = []
    DIH_ATOMS = []
    
    zmat_l = -1
    zmat_i = -1
    zmat_j = -1
    zmat_k = -1
    
    print ' finding dih_id and values from variables section of zmatrix '
    for line in iter(zmatrix.splitlines()) : 
	colvar = line.split('=')
        if( len( colvar) > 1 ):
            var_id = colvar[0].strip()
            if( var_id[0] == "D"  ):
                DIH_ID.append( colvar[0].strip() )
                DIH_VAL.append( float(colvar[1]) )
                DIH_ATOMS.append( [ zmat_l , zmat_i , zmat_j, zmat_k ] )
                
    print ' finding atoms of dih from zmatrix '
    line_n = 0
    for line in iter(zmatrix.splitlines()) : 
	col = line.split(',')
        line_n += 1
        
	if ( len(col) > 6 ):
            #print col[1],col[3]
            zmat_l = line_n
            zmat_i = int( col[1])
            zmat_j = int( col[3])
            zmat_k = int( col[5])
	    for dih_indx in  range(len(DIH_ID)):
		if(  col[6].strip() == DIH_ID[dih_indx].strip() ):
		    DIH_ATOMS[dih_indx] =   [ zmat_l  , zmat_i  , zmat_j , zmat_k ] 
                    
    return ( DIH_ID, DIH_VAL, DIH_ATOMS)
    


def com_zmatrix(com_name):
    """
    Get zmatrix from gaussian input file
    """

    # get lines 
    f = open(com_name,'r')
    zmatrix_lines = f.readlines()
    f.close()
    
    
    zmatrix = ""
    switch = -1 
    for line in zmatrix_lines: 
        col_space = line.split()
	colvar = line.split('=')        
        col = line.split(',')
        if ( switch == 1 and len(col[0]) < 2 ): 
            switch = 0
        if (  switch == 1 ):
	    zmatrix += line
            
        if (   len(col_space) == 2  ):
            if( col_space[0].strip() == "0" and col_space[1].strip() == "1" ):
            #if( isinstance(col[0], int) and isinstance(col[0], int) ):
                switch = 1
    
    return zmatrix



def tag_dih(strucC, RING_CONNECT, DIH_ID, DIH_VAL, DIH_ATOMS ):
    """
    Find dihedrals between two rings and set them loop for a torsional PES
    """

    verbose = True 
    
    DIH_TAG = []
    if ( verbose ): print " Tagging dihs "
    
    # Initialize all dihedrals to relax 
    for dih_indx in  range(len(DIH_ID)):
        DIH_TAG.append("relax")

    
    for dih_indx in  range(len(DIH_ID)):

        zmat_k = DIH_ATOMS[dih_indx][0]
        zmat_i = DIH_ATOMS[dih_indx][1]
        zmat_j = DIH_ATOMS[dih_indx][2]
        zmat_l = DIH_ATOMS[dih_indx][3]

        rn_k = strucC.ptclC[zmat_k].tagsDict["ring"]
        rn_i = strucC.ptclC[zmat_i].tagsDict["ring"]
        rn_j = strucC.ptclC[zmat_j].tagsDict["ring"]
        rn_l = strucC.ptclC[zmat_l].tagsDict["ring"]

        # Make sure inter ring connect not improper 
        if (  rn_i != rn_j and rn_l != rn_k ):
            DIH_TAG[dih_indx] = "fix"
	    
	    if ( verbose ):
		print "  Fixing dih ",DIH_ID[dih_indx], DIH_VAL[dih_indx], DIH_ATOMS[dih_indx]
	    
            if( rn_i > 0 and  rn_j  > 0   ):
                DIH_TAG[dih_indx] = "loop"
    
		if ( verbose ):
		    print "  Looping dih ",DIH_ID[dih_indx], DIH_VAL[dih_indx], DIH_ATOMS[dih_indx]
                    
    return DIH_TAG


def read_dihlist(dlist_name):
    
    
    DIH_ID = []
    DIH_VAL = []
    DIH_ATOMS = []
    DIH_TAG = []
    
    # get lines 
    f = open(dlist_name,'r')
    Lines = f.readlines()
    f.close()

    for line in Lines: 
        col = line.split()
        if( col[0] == "dih" and len(col) >= 11 ):
            dih_indx = col[1]
            DIH_TAG.append( col[2] )
            DIH_ID.append( col[3] )
            DIH_VAL.append(  float( col[4]) )
            a_l = int( col[5] ) - 1
            a_i = int( col[6] ) - 1
            a_j = int( col[7] ) - 1
            a_k = int( col[8] ) - 1
            DIH_ATOMS.append( [ a_l,a_i,a_j,a_k ] )
        
  
    return (  DIH_ID, DIH_VAL, DIH_TAG, DIH_ATOMS )
    
def write_dihlist(dlist_name,strucC, DIH_ID, DIH_VAL, DIH_TAG, DIH_ATOMS ):
    """
    Write dihidral list
    """
    verbose = True
    
    flist = open(dlist_name,"w")
    
    for dih_indx in  range(len(DIH_ID)):
        zmat_k = DIH_ATOMS[dih_indx][0]
        zmat_i = DIH_ATOMS[dih_indx][1]
        zmat_j = DIH_ATOMS[dih_indx][2]
        zmat_l = DIH_ATOMS[dih_indx][3]

        rn_k = strucC.ptclC[zmat_k].tagsDict["ring"]
        rn_i = strucC.ptclC[zmat_i].tagsDict["ring"]
        rn_j = strucC.ptclC[zmat_j].tagsDict["ring"]
        rn_l = strucC.ptclC[zmat_l].tagsDict["ring"]
        line_out = " dih %d %s %s %f %d %d %d %d %d %d %d %d \n"  % (dih_indx,DIH_TAG[dih_indx],DIH_ID[dih_indx], DIH_VAL[dih_indx],zmat_k,zmat_i,zmat_j,zmat_l ,rn_k,rn_i,rn_j,rn_l )
        if( verbose ): print line_out
        flist.write(line_out)
        
    flist.close()


def constrain_dih_zmatrix(zmatrix,DIH_ID,DIH_TAG,DIH_VAL,cent_angle,cent_indx):
    """
    Modify zmatrix to have a Constants dihedral value 
    """
    debug = False

    if( debug):
        print "zmatrix ",zmatrix
        
    con_zmatrix = ""
    # Prin all non constranied elments of zmatrix 
    for line in iter(zmatrix.splitlines()) :
        print_z = 1
	colvar = line.split('=')
        var_id = colvar[0].strip()
        for dih_indx in  range(len(DIH_ID)):
            if( DIH_ID[dih_indx] == var_id ):
                if(  DIH_TAG[dih_indx] != "relax"  ):
                    print_z = 0
                    break
            
        if( print_z ):
            con_zmatrix += "%s \n"% line 
            
            
    con_zmatrix += ' Constants: \n'
    # Print dihedral id's and angles as fixed
    for  dih_indx in  range(len(DIH_ID)):
	if (  DIH_ID[dih_indx] !=  DIH_ID[cent_indx]):
            if(  DIH_TAG[dih_indx] != "relax"  ):
                print dih_indx,DIH_ID[dih_indx],DIH_VAL[dih_indx]
                con_zmatrix +=   ' %s %8.4f   F  \n'%(DIH_ID[dih_indx],DIH_VAL[dih_indx])
	else:
	    con_zmatrix +=  '%s %8.4f  F  \n'%(DIH_ID[dih_indx],cent_angle)
	    
		
    con_zmatrix += ' \n'
    con_zmatrix += ' \n'

    if( debug):
        print "con_zmatrix  ",con_zmatrix
        sys.exit("debug con_zmatrix 1")
        

    return con_zmatrix

        
