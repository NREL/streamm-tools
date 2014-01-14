#! /usr/bin/env python
# IO for PDB file s

# Dr. Travis Kemper
# NREL
# Initial Date 05/07/2013
# travis.kemper@nrel.gov


#  Lists of atomic information for each atom in a structure  
#    ASYMB - list of atomic symbols
#    ELN   - list of atomic numbers
#    AMASS - list of atomic masses 

def atoms(Lines,options):
    import sys

    # declare arrays 
    ATYPE = []
    GTYPE = []
    RESID = []
    RESN = []
    CHARN = []
    R = []
    ASYMB = []
    EXCLUSIONS = []
    PDB_IND = []
    ATOM_NUMB = []
    # Read atoms section of pdb file
    
    debug = 0
    NA = 0

    ex_ids = options.exclusions.split()
    for indx_ex in range( len(ex_ids)):
        print " excluding :",ex_ids[indx_ex]
    
    for line in Lines :
        if ( line[0:4] == 'ATOM' ):
            atom_numb = int(line[5:11])
    print " allocating PDB_IND ", atom_numb
    for indx in range( atom_numb ):
        PDB_IND.append(atom_numb)
        
    PDB_IND.append(atom_numb)

    for line in Lines :
        if ( line[0:4] == 'ATOM' ):
            atom_numb = int(line[5:11])

            ATYPE_loc =  line[12:16] 
            GTYPE_loc = line[12:16] 
            ASYMB_loc =  line[76:78] 
            
            include_atom  = 1            
            for indx_ex in range( len(ex_ids)):
                if( ATYPE_loc == ex_ids[indx_ex] ):
                    include_atom  = 0
                if( GTYPE_loc == ex_ids[indx_ex] ):
                    include_atom  = 0
                if( ASYMB_loc == ex_ids[indx_ex] ):
                    include_atom  = 0
                    
            if( include_atom ):
                NA = NA + 1
                PDB_IND[atom_numb  ] = NA - 1
                ATYPE.append( ATYPE_loc )
                GTYPE.append( GTYPE_loc )
                RESID.append( line[17:20] )
                RESN.append( int(line[22:26]) )
                CHARN.append( int( line[22:26]) )
                R.append( [ float(line[30:38] ), float(line[38:46]), float(line[46:54]) ])
                ASYMB.append( line[76:78] )
                # CHARGES.append( line[78:80] )
                #CHARGES.append( float(0.0) )

                if( debug):
                    i = NA - 1
                    #print NA, ATYPE[i], ASYMB[i],RESN[i] ,  RESID[i] , CHARN[i] #, CHARGES[i]
                    print NA, i,atom_numb, ATYPE[i], ASYMB[i]
            else:
                ex_i =  int(line[5:11]) 
                EXCLUSIONS.append(ex_i)
                if( debug):
                    print " excluded ",ex_i,ATYPE_loc,GTYPE_loc,ASYMB_loc
                    
    debug = 0
    if( debug):
        for indx in range( len(PDB_IND)):
            pdb_i = PDB_IND[indx]
            print indx, pdb_i , ATYPE[pdb_i]
            
        sys.exit('atoms_pdb')

    return ( NA,ATYPE,GTYPE,RESID,RESN,CHARN,R,ASYMB,EXCLUSIONS ,PDB_IND ) 

def nablist(Lines,EXCLUSIONS,PDB_IND,ATYPE):
    import sys, numpy
    #
    debug = 0
    #
    maxnnab = 0
    numb_atoms = 0        
    for line in Lines :
        if ( line[0:6] == 'CONECT' ):
            maxnnab = maxnnab + 1
        if ( line[0:4] == 'ATOM' ):
            numb_atoms = int(line[5:11])
    maxnnab = maxnnab + 1
    maxnnab = 12*maxnnab 

    NBLIST = numpy.empty( maxnnab,  dtype=int )
    NBINDEX = numpy.empty( maxnnab,  dtype=int )


    if(debug):
        print " max nnab ",maxnnab

    NNAB = 0
    for line in Lines :
        if ( line[0:6] == 'CONECT' ):
            pdb_i = int(line[6:11]) # indexed from 0  -> NA -1
            include_i = 1
            for indx_ex in range( len(EXCLUSIONS)):
                ex_atom = int(EXCLUSIONS[indx_ex])
                if ( pdb_i == ex_atom ):
                    include_i = 0
                    if(debug): print " excluding atom  ", pdb_i, ex_atom
            if( include_i ):
                i = int(PDB_IND[pdb_i] )
                if(debug):
                    print " atom ",i+1,ATYPE[i],pdb_i,line[11:].split()
                ##    print " NNAB i ",NNAB,i
                NBINDEX[i] = NNAB + 1
                for line_j in line[11:].split() :
                    j = int(line_j.strip())
                    include_conect = 1
                    for indx_ex in range( len(EXCLUSIONS)):
                        ex_atom = int(EXCLUSIONS[indx_ex])
                        if ( j == ex_atom ):
                            include_conect = 0
                            if(debug): print " excluding connection ", j, ex_atom
                    if( include_conect ):
                        NNAB = NNAB + 1
                        if(debug): print "  for pdb index ", j
                        new_j = PDB_IND[j]
                        if(debug): print "    switching atomic index from  ", j," to ", new_j,ATYPE[new_j]  
                        NBLIST[NNAB] = new_j            
                        #if(debug):
                        #    print "NNAB  int(j) - 1 ",NNAB,int(j) - 1

    # Account for final atom position
    NBINDEX[i+1] =  NNAB + 1

    if(debug):
        sys.exit(" ex_atom debug in nablistd ")

    return (NBLIST, NBINDEX)
