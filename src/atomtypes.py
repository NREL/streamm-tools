#! /usr/bin/env python
"""
Asign atom types for force-field topologies 
"""

# Dr. Travis Kemper
# NREL
# Initial Date 05/07/2013
# travis.kemper@nrel.gov

import sys
import topology as top 

def oplsaa( update_chr ,strucC , ring_nblist, ring_nbindex,cov_nblist, cov_nbindx ):
    """
    Set OPLSaa atom types 
    """
    
    
    for pid_i, ptclObj_i  in strucC.ptclC:
        NNAB = top.calc_nnab(pid_i,cov_nbindx)
        ELCNT = top.calc_elcnt(pid_i,strucC,cov_nblist,cov_nbindx)
        if ptclObj_i.tagsDict["number"] == 6 :
           
            # simple guess based on coordination 
            if int(NNAB) == 4 :
                ptclObj_i.tagsDict["fftype"] = 'CT' # Alkane
                # refine guess based on nieghbors 
                if( ELCNT[1] == 4 ):                              # Methane 
                    ptclObj_i.tagsDict["fftype"] = 'CT'
                    if( update_chr ): ptclObj_i.charge =  -0.24
                elif( ELCNT[6] == 1 and ELCNT[1] == 3):         # Methyl 
                    ptclObj_i.tagsDict["fftype"] = 'CT'
                    if( update_chr ): ptclObj_i.charge =  -0.18
                elif(   ELCNT[6] == 2 and ELCNT[1] == 2 ):      
                    ptclObj_i.tagsDict["fftype"] = 'CT'
                    if( update_chr ): ptclObj_i.charge =  -0.12
                elif( ELCNT[6] == 3 and ELCNT[1] == 1 ):         # Alkane 
                    ptclObj_i.tagsDict["fftype"] = 'CT'
                    if( update_chr ): ptclObj_i.charge =  -0.06
                elif(  ELCNT[6] == 4 ):
                    ptclObj_i.tagsDict["fftype"] = 'CT'
                    if( update_chr ): ptclObj_i.charge =  0.0
                elif(  ELCNT[1] == 1 and ELCNT[8] == 2 ):       #acetal
                    ptclObj_i.tagsDict["fftype"] = 'CO'
                    if( update_chr ): ptclObj_i.charge =  -0.4

            if  int(NNAB) == 3 :
                r_numb = ptclObj_i.tagsDict["ring"]
                if( r_numb != 0 ):
                    nring = top.calc_nnab(r_numb,ring_nbindex)
                    nring_mod = nring % 2     # to tell 6 adn 5 member rings apart
                    ptclObj_i.tagsDict["fftype"] = 'CA' # Aromatic C
                    if( ELCNT[6] == 2 and ELCNT[1] == 1 ):           # "Aromatic C"  
                        ptclObj_i.tagsDict["fftype"] = 'CA'
                        if( update_chr ): ptclObj_i.charge =  -0.1150
                    if( ELCNT[6] == 3  and nring_mod == 0 ):            # "Naphthalene Fusion C"
                        ptclObj_i.tagsDict["fftype"] = 'CA'
                        if( update_chr ): ptclObj_i.charge =  0.0
                    if( nring == 5 or nring == 8 or nring == 12 ):   # Ring with 5 membered parts
                        if( ELCNT[6] == 3 ):            # Fused 5-6 rings 
                            ptclObj_i.tagsDict["fftype"] = 'CB'
                            if( update_chr ): ptclObj_i.charge =  0.0
                        if( ELCNT[7] == 1 or ELCNT[16] == 1  ):      # Pyrrole or thiophene 
                            if( ELCNT[6] == 1 and ELCNT[1] == 1 ):             #imidazole C2
                                ptclObj_i.tagsDict["fftype"] = 'CW'   
                                if( update_chr ): ptclObj_i.charge =  0.22
                            if( ELCNT[6] == 1 and ELCNT[8] == 1 ):             # Not sure about this one
                                ptclObj_i.tagsDict["fftype"] = 'CW'   
                                if( update_chr ): ptclObj_i.charge =  0.22
                        if(  ELCNT[7] == 2 ):             #imidazole C2
                            ptclObj_i.tagsDict["fftype"] = 'CR'   
                            if( update_chr ): ptclObj_i.charge =  0.22
                        if(  ELCNT[7] == 1 and ELCNT[16] == 1 ):             #imidazole C2
                            ptclObj_i.tagsDict["fftype"] = 'CR'
                            if( update_chr ): ptclObj_i.charge =  0.22

                    if( ELCNT[6] == 2 and  ELCNT[8] == 1  ):            # CTD strangeness ... 
                        ptclObj_i.tagsDict["fftype"] = 'CW'
                        if( update_chr ): ptclObj_i.charge =  0.0
                        
                    debug = 0
                    if( debug  ):
                        print " r_numb ",nring,r_numb, nring_mod, ELCNT[6], ELCNT[7], ptclObj_i.tagsDict["fftype"]
                        # sys.exit(' debug ')
                    
                else: # not aromatic
                    if( ELCNT[6] == 2 and ELCNT[1] == 1 ):          # diene 
                        ptclObj_i.tagsDict["fftype"] = 'C='
                        if( update_chr ): ptclObj_i.charge =  0.0
                    elif( ELCNT[6] == 2 and ELCNT[8] == 1 ):          # Benzophenone
                        ptclObj_i.tagsDict["fftype"] = 'C'
                        if( update_chr ): ptclObj_i.charge =  0.7
                    elif( ELCNT[6] == 1 and ELCNT[8] == 2 ):
                        ptclObj_i.tagsDict["fftype"] = 'C' # pmma
                        if( update_chr ): ptclObj_i.charge =  0.7

            if int(NNAB) == 2 :
                ptclObj_i.tagsDict["fftype"] = 'C:'    # Allene
                if( ELCNT[6] == 1 and ELCNT[7] == 1 ):   # "Benzonitrile -CN"  
                    ptclObj_i.tagsDict["fftype"] = 'CZ'
                    
            if int(NNAB) == 1 :
                ptclObj_i.tagsDict["fftype"] = '' # Aromatic C
                print " WARNING!!! carbon index ",atom_i," bonded to single atom "
        #
        # label oxygens
        #
        if( ptclObj_i.tagsDict["number"] == 8 ):
            if int(NNAB) == 1 :
                ptclObj_i.tagsDict["fftype"] = 'O' # double bonded
                if( update_chr ): ptclObj_i.charge =  -0.5
            if int(NNAB) == 2 :
                ptclObj_i.tagsDict["fftype"] = 'OS' # ether
                if( update_chr ): ptclObj_i.charge =  -0.5
                if( ELCNT[1] == 1 ):
                    ptclObj_i.tagsDict["fftype"] = 'OH' # Alcohol
                    if( update_chr ): ptclObj_i.charge =   -0.6830
                if( ELCNT[8] == 1 ):
                    ptclObj_i.tagsDict["fftype"] = 'O2' # Carboxylate
                    if( update_chr ): ptclObj_i.charge =   -0.800
                if( ELCNT[16] == 1 ):
                    ptclObj_i.tagsDict["fftype"] = 'OY' # Sulfoxide
                    if( update_chr ): ptclObj_i.charge =   -0.4200
            if int(NNAB) == 3 :
                if( ELCNT[7] == 1  ):
                    ptclObj_i.tagsDict["fftype"] = 'ON'
                    if( update_chr ): ptclObj_i.charge =  -0.118
                        
        #
        # label nitrogens 
        #
        if ptclObj_i.tagsDict["number"] == 7 :
            if int(NNAB) == 3 :      # amide
                ptclObj_i.tagsDict["fftype"] = 'N' 
                if( ELCNT[1] == 3 ): # Ammonia NH3"   
                    ptclObj_i.tagsDict["fftype"] = 'NT'
                if( ELCNT[1] == 2 ): # -NH2
                    ptclObj_i.tagsDict["fftype"] = 'N2'
                if( ELCNT[1] == 1 ): # -NH2
                    ptclObj_i.tagsDict["fftype"] = 'N2'
                    
            r_numb = ptclObj_i.tagsDict["ring"]
            if( r_numb != 0 ):
                nring = top.calc_nnab(r_numb,ring_nbindex)
                nring_mod = nring % 2     # to tell 6 adn 5 member rings apart
                if ( nring == 5 or nring == 8 or nring == 12 ):   # Ring with 5 membered parts
                    if ( int(NNAB) == 2 ):
                        ptclObj_i.tagsDict["fftype"] = 'NC'      # "Imidazole N3"     
                    if  int(NNAB) == 3 :
                        ptclObj_i.tagsDict["fftype"] = 'NA'      # "Imidazole N1"      
                else:
                    if ( int(NNAB) == 2  ):
                        ptclObj_i.tagsDict["fftype"] = 'NB'
            if ( int(NNAB) == 1 and ELCNT[6] == 1 ) :      # Nitrile 
                ptclObj_i.tagsDict["fftype"] = 'NZ'
                

        #
        # label sulfurs
        #
        if( ptclObj_i.tagsDict["number"] == 16 ):
            if int(NNAB) == 2 :
                ptclObj_i.tagsDict["fftype"] = 'S'   #  Thioether RSR (UA)
                if( update_chr ): ptclObj_i.charge =-0.4350
                if( ELCNT[1] == 1  ):
                    ptclObj_i.tagsDict["fftype"] = 'SH'
                    if( update_chr ): ptclObj_i.charge =  -0.335
                if( ELCNT[1] == 2  ):
                    ptclObj_i.tagsDict["fftype"] = 'SH'
                    if( update_chr ): ptclObj_i.charge =  -0.470
                        
            if ( int(NNAB) == 4 and ELCNT[8] == 2 and ELCNT[7] == 1 ):
                ptclObj_i.tagsDict["fftype"] = 'SY'   #  "Sulfonamide -SO2N<" 
            if ( int(NNAB) == 3 and ELCNT[8] == 1 ):
                ptclObj_i.tagsDict["fftype"] = 'SZ'   #  Sulfoxide   
                           
        # Label chlorine 
        if ( ptclObj_i.tagsDict["number"] == 17):
            ptclObj_i.tagsDict["fftype"] = 'Cl'
            if( update_chr ): ptclObj_i.charge =  -0.2

        # Label Fluorine   
        if ( ptclObj_i.tagsDict["number"] == 9):
            ptclObj_i.tagsDict["fftype"] = 'F'
            if( update_chr ): ptclObj_i.charge =  -0.2057

    #
    # label hydrogens
    #

    for pid_i, ptclObj_i  in strucC.ptclC:
        NNAB = top.calc_nnab(pid_i,cov_nbindx)
        ELCNT = top.calc_elcnt(pid_i,strucC,cov_nblist,cov_nbindx)
        if( ptclObj_i.tagsDict["number"] == 1 ):
            if ( NNAB > 1  ):
                sys.exit(' over coordinated H')
            pid_j = cov_nblist[ cov_nbindx[pid_i] ]
            ptclObj_j = strucC.ptclC[pid_j]
            ELCNT_j =  top.calc_elcnt(pid_j,strucC,cov_nblist,cov_nbindx)
            if ( ptclObj_j.tagsDict["fftype"]== 'CA' ):
                ptclObj_i.tagsDict["fftype"] = 'HA' #
                if( update_chr ): ptclObj_i.charge =  0.115
            if ( ptclObj_j.tagsDict["fftype"]== 'CT' ):
                ptclObj_i.tagsDict["fftype"] = 'HC' #
                if( update_chr ): ptclObj_i.charge =  0.06
            if ( ptclObj_j.tagsDict["fftype"]== 'CW' ):
                ptclObj_i.tagsDict["fftype"] = 'HA' #
            if ( ptclObj_j.tagsDict["fftype"]== 'CS' ):
                ptclObj_i.tagsDict["fftype"] = 'HA' #
            if( ELCNT_j[6] == 2 and ELCNT_j[8] == 1  ): # "Ester -OCH<"
                ptclObj_i.tagsDict["fftype"] = 'H1'
                if( update_chr ): ptclObj_j.charge =  0.03
            if ( ptclObj_j.tagsDict["fftype"]== 'SH' ):
                ptclObj_i.tagsDict["fftype"] = 'HS' #
                if( update_chr ): ptclObj_i.charge = 0.1550
            if ( ptclObj_j.tagsDict["number"] == 7 ):
                ptclObj_i.tagsDict["fftype"] = 'H' #
            if ( ptclObj_j.tagsDict["fftype"]== 'NA' ):
                ptclObj_i.tagsDict["fftype"] = 'H' #
            if ( ptclObj_j.tagsDict["fftype"]== 'N2' ):
                ptclObj_i.tagsDict["fftype"] = 'H' #
            if ( ptclObj_j.tagsDict["fftype"]== 'OH' ):
                ptclObj_i.tagsDict["fftype"] = 'HO' #
                    
                    
    # relabel based on neighbors
    for pid_i, ptclObj_i  in strucC.ptclC:
        NNAB = top.calc_nnab(pid_i,cov_nbindx)
        ELCNT = top.calc_elcnt(pid_i,strucC,cov_nblist,cov_nbindx)
        N_o = cov_nbindx[ pid_i ]
        N_f = cov_nbindx[ pid_i + 1 ] - 1
        nring_i = ptclObj_i.tagsDict["ring"]
        if ( ptclObj_i.tagsDict["fftype"] == 'CA' ):
            for indx_j in range( N_o,N_f+1):
                pid_j = cov_nblist[indx_j]
                ptclObj_j = strucC.ptclC[pid_j]
                ELCNT_j =  top.calc_elcnt(pid_j,strucC,cov_nblist,cov_nbindx)
                nring_j = ptclObj_j.tagsDict["ring"]
                if ( ptclObj_j.tagsDict["fftype"]== 'CW' ):
                    if( ELCNT[1] == 1 and ELCNT[6] == 2 ):
                        ptclObj_i.tagsDict["fftype"] = 'CS'
                if(  ELCNT[6] == 2 and nring_i == nring_j and ptclObj_j.tagsDict["number"] != 6 ):          # fussed 
                    if( debug):
                        print ' CB ',pid_i+1, ELCNT[6] ,  nring_i , nring_j 
                    ptclObj_i.tagsDict["fftype"] = 'CB'
            if ( ptclObj_i.tagsDict["fftype"] == 'CS' ):
                for indx_j in range( N_o,N_f+1):
                    pid_j = cov_nblist[indx_j]
                    ptclObj_j = strucC.ptclC[pid_j]
                    ELCNT_j = top.calc_elcnt(pid_j,strucC,cov_nblist,cov_nbindx)
                    if (  ptclObj_j.tagsDict["fftype"]== 'CA' ):
                        if ( ELCNT_j[6] == 3 ):
                            ptclObj_j.tagsDict["fftype"]= 'CS'
                        if ( ELCNT_j[6] == 2 and ELCNT_j[1] == 1 ):
                            ptclObj_j.tagsDict["fftype"]= 'CS'

        if ptclObj_i.tagsDict["fftype"] == 'CT' :
            for indx_j in range( N_o,N_f+1):
                pid_j = cov_nblist[indx_j]
                ptclObj_j = strucC.ptclC[pid_j]
                if ( ptclObj_j.tagsDict["fftype"]== 'O' ):
                    if( ELCNT[1] == 3 and ELCNT[8] == 1 ):
                        ptclObj_i.tagsDict["fftype"] = 'C3' # Methyloxide
                        
                        
    # Find Ring linkers
    debug = 0
    for pid_i, ptclObj_i  in strucC.ptclC:
        NNAB = top.calc_nnab(pid_i,cov_nbindx)
        ELCNT = top.calc_elcnt(pid_i,strucC,cov_nblist,cov_nbindx)
        N_o = cov_nbindx[ pid_i ]
        N_f = cov_nbindx[ pid_i + 1 ] - 1
        if ( ptclObj_i.tagsDict["fftype"] == 'CA' or  ptclObj_i.tagsDict["fftype"] == 'CB'  or  ptclObj_i.tagsDict["fftype"] == 'CW' ) :
            if( NNAB == 3 ):  # if sp2
                for indx_j in range( N_o,N_f+1):
                    pid_j = cov_nblist[indx_j]
                    ptclObj_j = strucC.ptclC[pid_j]
                    if ( ptclObj_j.tagsDict["fftype"]== 'CA' or  ptclObj_j.tagsDict["fftype"]== 'CB'  or  ptclObj_j.tagsDict["fftype"]== 'CW' )  :
                        NNAB_j = top.calc_nnab(pid_j,cov_nbindx)
                        if( NNAB_j == 3 ):  # if sp2
                            if(  ptclObj_i.tagsDict["ring"] !=  ptclObj_j.tagsDict["ring"] ): # ring linker
                                if(debug):
                                    print pid_i,' and ',pid_j,' link ',ptclObj_i.tagsDict["fftype"],ptclObj_i.tagsDict["fftype"]
                                    print ' ', ptclObj_i.tagsDict["ring"] , ptclObj_j.tagsDict["ring"]
                                ptclObj_i.tagsDict["fftype"] = 'C!'
                                ptclObj_j.tagsDict["fftype"]= 'C!'                                

    # relabel based on neighbors
    debug = 0
    for pid_i, ptclObj_i  in strucC.ptclC:
        if ( ptclObj_i.tagsDict["fftype"] == 'C!' ):
            r_numb = ptclObj_i.tagsDict["ring"]
            nring =  top.calc_nnab(r_numb,ring_nbindex)
            ELCNT = top.calc_elcnt(pid_i,strucC,cov_nblist,cov_nbindx)
            if( debug ):
                print ' linker found ',atom_i
                print ' nring,ELCNT[6],ELCNT ',nring,ELCNT[6],ELCNT 
            if( nring == 5 or nring == 8 or nring == 12 ):   # Ring with 5 membered parts
                if( ELCNT[7] == 1 or ELCNT[16] == 1  ):      # Pyrrole or thiophene 
                    if( ELCNT[6] == 2  ):             #imidazole C2
                        N_o = cov_nbindx[ pid_i ]
                        N_f = cov_nbindx[ pid_i + 1 ] - 1
                        for indx_j in range( N_o,N_f+1):
                            pid_j = cov_nblist[indx_j]
                            ptclObj_j =  strucC.ptclC[pid_j]
                            if( ptclObj_j.tagsDict["fftype"]== 'CA' ):
                                if( debug ): print " changing ",ptclObj_j.tagsDict["fftype"], " to CS"
                                ptclObj_j.tagsDict["fftype"]= 'CS'   
                                if( update_chr ): ptclObj_i.charge =  0.22

    if(debug):   sys.exit('linkers')
  
    debug = 0
    if(debug):
        for pid_i, ptclObj_i  in strucC.ptclC:
            print pid_i,ptclObj_i.tagsDict["fftype"],ptclObj_i.tagsDict["number"]
        sys.exit('debug')

    # Check for unidentified atoms
    for pid_i, ptclObj_i  in strucC.ptclC:
        if ( ptclObj_i.tagsDict["fftype"] == '?' ):
            print ' atom ', pid_i , ptclObj_i.tagsDict["number"],' unknow '
            sys.exit(' Unknow atom ')
         
    return (strucC)

def set_ptmatypes( update_chr,strucC , ring_nblist, ring_nbindex,cov_nblist, cov_nbindx ):
    """
    Set nitroxyl atom types for TEMPO 
    """
    

    set_resnames = False 
    resname_i = 'PTMA'
    methyl_C = "C6"

    debug = 0
     
    
    # Set default charges, which will be over writen in some cases 
    for pid_i, ptclObj_i  in strucC.ptclC:

        if(debug): print ptclObj_i.tagsDict["symbol"]
        #if ptclObj_i.tagsDict["fftype"] == 'HC' :
        #    if( update_chr ): ptclObj_i.charge =  0.06
        #if ptclObj_i.tagsDict["fftype"] == 'CT' :
        #    if( update_chr ): ptclObj_i.charge =  -0.18
        if ptclObj_i.tagsDict["symbol"].strip() == 'LP' :
            ptclObj_i.tagsDict["fftype"] = 'LP'
            if(set_resnames): ptclObj_i.tagsDict["resname"] = resname_i
            ptclObj_i.mass = 9.0 
            if( update_chr ): ptclObj_i.charge =  -0.11
        
    for pid_i, ptclObj_i  in strucC.ptclC:
        NNAB = top.calc_nnab(pid_i,cov_nbindx)
        ELCNT = top.calc_elcnt(pid_i,strucC,cov_nblist,cov_nbindx)
        #
        # label nitrogens 
        #
        if ( ptclObj_i.tagsDict["number"] == 7 ):
            if( ELCNT[8] == 1 and ELCNT[6] == 2  ):
                Ni_o = cov_nbindx[ pid_i ]
                Ni_f = cov_nbindx[ pid_i+1  ] - 1
                nitroxide = 0
                cation = 0
                for indx_i in range( Ni_o,Ni_f+1):
                    pid_j = cov_nblist[indx_i]
                    ptclObj_j = strucC.ptclC[pid_j]
                    if ( ptclObj_j.tagsDict["number"] == 8 ):
                        NNAB_j = top.calc_nnab(pid_j,cov_nbindx)
                        ELCNT_j = top.calc_elcnt(pid_j,strucC,cov_nblist,cov_nbindx)
                        if( NNAB_j == 1 ):
                            nitroxide = 1
                            cation = 1
                        elif( NNAB_j == 3 and  ELCNT_j[0] == 2 ):
                            nitroxide = 1
                            cation = 0
                        if(debug): print ' oxygen ', NNAB_j,ELCNT_j[0]
                if( nitroxide ):
                    if(debug): print " nitroxide found ",cation
                    
                    if( cation == 0 ):
                        residue = 'PTMA'
                        q_NN    =  0.076
                        q_ON    = -0.209
                        q_CNN   = 0.558
                        q_CTN   = -0.289
                        q_HN    = 0.06
                        q_CT2   = -0.404
                        q_HCT2  = 0.06
                        # same as dist 1
                        q_Cp    = 0.489
                        q_H1    = 0.05
                        q_C     = 0.679
                        q_OS    = -0.503
                        q_O     = -0.503
                        q_Cb1   = 0.159
                        q_CH2   = -0.120
                        q_CH3   = -0.310
                        q_H3    = 0.06
                    elif( cation == 1 ):
                        residue = 'PTMc'
                        # Use cation charges
                        q_NN    =  0.635
                        q_ON    = -0.195
                        q_CNN   = 0.180
                        q_CTN   = -0.143
                        q_HN    = 0.06
                        q_CT2   = -0.235
                        q_HCT2  = 0.06
                        q_Cp    = 0.435
                        q_H1    = 0.065
                        q_C     = 0.757
                        q_OS    = -0.451
                        q_O     = -0.574
                        q_Cb1   = 0.126
                        q_CH2   = -0.120
                        q_CH3   = -0.256
                        q_H3    = 0.06

                    ptclObj_i.tagsDict["fftype"] = 'NN'
                    if( update_chr ): ptclObj_i.charge = q_NN
                    if(set_resnames): ptclObj_i.tagsDict["resname"] = resname_i
                    for indx_i in range( Ni_o,Ni_f+1):
                        pid_j = cov_nblist[indx_i]
                        ptclObj_j = strucC.ptclC[pid_j]
                        if ( ptclObj_j.tagsDict["number"] == 8 ):
                            ptclObj_j.tagsDict["fftype"] = 'ON'
                            if( update_chr ): ptclObj_j.charge = q_ON
                                    
                            if(set_resnames): ptclObj_j.tagsDict["resname"] = resname_i
                        if ( ptclObj_j.tagsDict["number"] == 6 ):
                            ptclObj_j.tagsDict["fftype"]= 'CT'
                            if( update_chr ): ptclObj_j.charge = q_CNN
                            if(set_resnames): ptclObj_j.tagsDict["resname"] = resname_i
                            Nj_o = cov_nbindx[ pid_j ]
                            Nj_f = cov_nbindx[ pid_j +1  ] - 1
                            for indx_j in range(Nj_o,Nj_f+1):
                                pid_k = cov_nblist[indx_j]
                                ptclObj_k = strucC.ptclC[pid_k]
                                if(set_resnames):  ptclObj_k.tagsDict["resname"]  = resname_i
                                ELCNT_k = top.calc_elcnt(pid_k,strucC,cov_nblist,cov_nbindx)
                                if( ptclObj_k.tagsDict["number"] == 6 ):
                                    if ( ELCNT_k[1] == 3 ):
                                        if( update_chr ):  ptclObj_k.charge =  q_CTN
                                        if(set_resnames):  ptclObj_k.tagsDict["resname"]  = resname_i
                                        Nk_o = cov_nbindx[ pid_k ]
                                        Nk_f = cov_nbindx[ pid_k+1  ] - 1
                                        for indx_k in range(Nk_o,Nk_f+1):
                                            pid_l = cov_nblist[indx_k]
                                            ptclObj_l = strucC.ptclC[pid_l]
                                            if( ptclObj_l.tagsDict["number"] == 1):
                                                ptclObj_l.tagsDict["fftype"] = 'HN'
                                                if( update_chr ): ptclObj_l.charge=  q_HN
                                                if(set_resnames):  ptclObj_l.tagsDict["3"] = resname_i
                                            
                                    if ( ELCNT_k[1] == 2  and ELCNT_k[1] == 2 ):       
                                        if( update_chr ): ptclObj_k.charge =  q_CT2
                                        if(set_resnames): ptclObj_k.tagsDict["resname"] = resname_i
                                        Nk_o = cov_nbindx[ pid_k ]
                                        Nk_f = cov_nbindx[ pid_k+1  ] - 1
                                        for indx_k in range(Nk_o,Nk_f+1):
                                            pid_l = cov_nblist[indx_k]
                                            ptclObj_l = strucC.ptclC[pid_l]
                                            if(set_resnames):  ptclObj_l.tagsDict["resname"] = resname_i
                                            if( ptclObj_l.tagsDict["number"] == 1 ):
                                                if( update_chr ): ptclObj_l.charge =  q_HCT2
                                                if(set_resnames): ptclObj_l.tagsDict["resname"] = resname_i
                                            if( ptclObj_l.tagsDict["number"] == 6 and pid_l != pid_k and pid_l != pid_j ):
                                                ptclObj_l.tagsDict["fftype"] = 'CT'
                                                if( update_chr ): ptclObj_l.charge =  q_Cp
                                                if(set_resnames):  ptclObj_l.tagsDict["resname"] = resname_i
                                                Nl_o = cov_nbindx[ pid_l ]
                                                Nl_f = cov_nbindx[ pid_l+1  ] - 1
                                                for indx_l in range(Nl_o,Nl_f+1):
                                                    pid_m = cov_nblist[indx_l]
                                                    ptclObj_m = strucC.ptclC[pid_m]
                                                    if(set_resnames): ptclObj_m.tagsDict["resname"]  = resname_i
                                                    if( ptclObj_m.tagsDict["number"]  == 1 ):
                                                        ptclObj_m.tagsDict["fftype"] = 'H1'
                                                        if( update_chr ): ptclObj_m.charge =  q_H1
                                                        if(set_resnames): ptclObj_m.tagsDict["resname"]  = resname_i
                                                    if( ptclObj_m.tagsDict["number"]  == 8 ):
                                                        ptclObj_m.tagsDict["fftype"]  = 'OS'
                                                        if( update_chr ): ptclObj_m.charge  =  q_OS
                                                        if(set_resnames): ptclObj_m.tagsDict["resname"] = resname_i
                                                        Nm_o = cov_nbindx[ pid_m ]
                                                        Nm_f = cov_nbindx[ pid_m + 1 ] - 1
                                                        for indx_m in range(Nm_o,Nm_f+1):
                                                            pid_n = cov_nblist[indx_m]
                                                            ptclObj_n = strucC.ptclC[pid_n]
                                                            if(set_resnames): ptclObj_n.tagsDict["resname"] = resname_i
                                                            if( ptclObj_n.tagsDict["number"] == 6 and pid_n != pid_l ):
                                                                ptclObj_n.tagsDict["fftype"]  = 'C'
                                                                if( update_chr ): ptclObj_n.charge   =  q_C
                                                                if(set_resnames): ptclObj_n.tagsDict["resname"]  = resname_i
                                                                Nn_o = cov_nbindx[ pid_n ]
                                                                Nn_f = cov_nbindx[ pid_n+1  ] - 1
                                                                for indx_n in range(Nn_o,Nn_f+1):
                                                                    pid_o = cov_nblist[indx_n]
                                                                    ptclObj_o = strucC.ptclC[pid_o]
                                                                    if( ptclObj_o.tagsDict["number"] == 8 and pid_o != pid_m ):
                                                                        ptclObj_o.tagsDict["fftype"] = 'O'
                                                                        if( update_chr ): ptclObj_o.charge  =  q_O
                                                                        if(set_resnames): ptclObj_o.tagsDict["resname"] = resname_i
                                                                    if( ptclObj_o.tagsDict["number"]  == 6 ):
                                                                        ptclObj_o.tagsDict["fftype"] = 'CT'
                                                                        if( update_chr ): ptclObj_o.charge=  q_Cb1
                                                                        if(set_resnames): ptclObj_o.tagsDict["resname"]= 'BCK' #resname_i
                                                                        No_o = cov_nbindx[ pid_o ]
                                                                        No_f = cov_nbindx[ pid_o+1  ] - 1
                                                                        for indx_o in range(No_o,No_f+1):
                                                                            pid_p = cov_nblist[indx_o]
                                                                            ptclObj_p = strucC.ptclC[pid_p]
                                                                            ELCNT_p = top.calc_elcnt(pid_p,strucC,cov_nblist,cov_nbindx)
                                                                            if( ptclObj_p.tagsDict["number"]  == 6 and pid_p != pid_n ):
                                                                                if( ELCNT_p[1] == 2 ):
                                                                                    ptclObj_p.tagsDict["fftype"] = 'CT'
                                                                                    if( update_chr ): ptclObj_p.charge =  q_CH2
                                                                                    if(set_resnames): ptclObj_p.tagsDict["resname"] = 'BCK' #resname_i
                                                                                    Np_o = cov_nbindx[ pid_p ]
                                                                                    Np_f = cov_nbindx[ pid_p+1  ] - 1
                                                                                    for indx_p in range(Np_o,Np_f+1):
                                                                                        pid_q = cov_nblist[indx_p]
                                                                                        ptclObj_q = strucC.ptclC[pid_q]
                                                                                        if( ptclObj_q.tagsDict["number"]  == 1 ):
                                                                                            #ATYPE[pid_q] = 'H1'
                                                                                            #if( update_chr ): CHARGES[pid_q] = q_H3
                                                                                            if(set_resnames): ptclObj_q.tagsDict["resname"] = 'BCK' # resname_i
                                                                                        
                                                                                if( ptclObj_p.tagsDict["gtype"].strip()  == "C36" or ptclObj_p.tagsDict["gtype"].strip()  == "C77"  ): #or  # ELCNT_l[1] == 3 ):
                                                                                    
                                                                                    if( update_chr ): ptclObj_p.charge = q_CH3
                                                                                    if(set_resnames): ptclObj_p.tagsDict["resname"] = resname_i
                                                                                    Np_o = cov_nbindx[pid_p]
                                                                                    Np_f = cov_nbindx[pid_p+1] - 1
                                                                                    for indx_p in range( Np_o,Np_f+1):
                                                                                        pid_q = cov_nblist[indx_p]
                                                                                        ptclObj_q = strucC.ptclC[pid_q]
                                                                                        if( ptclObj_q.tagsDict["number"]  == 1 ):
                                                                                            #ATYPE[pid_q] = 'H1'
                                                                                            if( update_chr ): ptclObj_q.charge = q_H3
                                                                                            if(set_resnames): ptclObj_q.tagsDict["resname"] = resname_i
                                                                                            
                                                    

    debug=False
    if( debug ):
        for pid_i, ptclObj_i  in strucC.ptclC:
            print pid_i ,ptclObj_i.tagsDict["number"] ,ptclObj_i.tagsDict["gtype"],ptclObj_i.tagsDict["fftype"],ptclObj_i.charge ,ptclObj_i.tagsDict["resname"]
        sys.exit('ptma_types ')
    #
    ## Set charge groups 
    #n_groups = 0
    #    
    ## set heavy atom groups
    #for i in range(NA):
    #    q_g = CHARN[pid_i]
    #    if [ q_g > n_groups ]: n_groups =  q_g 
    #if( options.verbose ):
    #    print " "
    ## set heavy atom groups
    #for i in range(NA):
    #    if( CG_SET[i] ):
    #        NNAB = calc_nnab(i,cov_nblist,cov_nbindx)
    #        ELCNT = calc_elcnt(i,ELN,cov_nblist,cov_nbindx)
    #        N_o = cov_nbindx[i]
    #        N_f = cov_nbindx[i+1] - 1
    #        if( ATYPE[i] == 'NN' ):  # Methyl 
    #            n_groups = CHARN[i]
    #            CG_SET[i] = 0
    #            CHARN[i] = n_groups
    #            for j_indx in range( N_o,N_f+1):
    #                j = cov_nblist[j_indx]
    #                #if ( ASYMB[j] == 'H' ):
    #                CHARN[j] = n_groups
    #                CG_SET[j] = 0
    #                    
                        
    return ( strucC )

def set_pmmatypes(update_chr,strucC , ring_nblist, ring_nbindex,cov_nblist, cov_nbindx ):
    """
    Set atom types for PMMA monomer
    """

    methyl_C = "C6"
    debug = False
    set_resnames = False

    # Version 2
    q_OMe = -0.251
    q_OHMe = 0.13
    q_C = 0.894
    q_O = -0.611
    q_OS = -0.313
    q_BC =  -0.012
    q_CH2 = -0.12
    q_CH3 = -0.46
    q_H3 = 0.121

    # Version 1
    q_OMe = -0.041
    q_OHMe = 0.06
    q_C = 0.894
    q_O = -0.611
    q_OS = -0.313
    q_BC =  -0.012
    q_CH2 = -0.12
    q_CH3 = -0.277
    q_H3 = 0.06

    residue = 'PMMA' 

    # Find pmma segments and label charges
    for pid_i, ptclObj_i  in strucC.ptclC:
        NNAB = top.calc_nnab(pid_i,cov_nbindx)
        ELCNT = top.calc_elcnt(pid_i,strucC,cov_nblist,cov_nbindx)
        N_o = cov_nbindx[pid_i]
        N_f = cov_nbindx[pid_i+1] - 1
        if ptclObj_i.tagsDict["fftype"] == 'OS' :
            pmma_oxygen = 0
            for indx in range( N_o,N_f+1):
                pid_j = cov_nblist[indx]
                ptclObj_j = strucC.ptclC[pid_j]
                NNAB_j =  top.calc_nnab(pid_j,cov_nbindx)
                ELCNT_j = top.calc_elcnt(pid_j,strucC,cov_nblist,cov_nbindx)
                if( ELCNT_j[1] == 3 ):
                    pmma_oxygen = 1
            if( pmma_oxygen == 1 ):
                if( update_chr ): ptclObj_i.charge = q_OS
                if(set_resnames): ptclObj_i.tagsDict["resname"] = residue
                for indx in range( N_o,N_f+1):
                    pid_j = cov_nblist[indx]
                    ptclObj_j = strucC.ptclC[pid_j]
                    if( ptclObj_j.tagsDict["fftype"] == 'CT' ):
                        if( update_chr ):
                            ptclObj_j.charge = q_OMe
                        if(set_resnames): ptclObj_j.tagsDict["resname"] = residue
                        Nj_o = cov_nbindx[pid_j]
                        Nj_f = cov_nbindx[pid_j+1] - 1
                        for indx_j in range( Nj_o,Nj_f+1):
                            pid_k = cov_nblist[indx_j]
                            ptclObj_k = strucC.ptclC[pid_k]
                            if( ptclObj_k.tagsDict["fftype"] == 'HC' ):
                                ptclObj_k.tagsDict["fftype"]  = 'H1'
                                if( update_chr ): ptclObj_k.charge = q_OHMe
                                if(set_resnames):  ptclObj_k.tagsDict["resname"] = residue
                                #
 
                    if( ptclObj_j.tagsDict["fftype"] == 'C' ):
                        if( update_chr ): ptclObj_j.charge = q_C
                        if(set_resnames):  ptclObj_j.tagsDict["resname"] = residue
                        Nj_o = cov_nbindx[pid_j]
                        Nj_f = cov_nbindx[pid_j+1] - 1
                        for indx_j in range( Nj_o,Nj_f+1):
                            pid_k = cov_nblist[indx_j]
                            ptclObj_k = strucC.ptclC[pid_k]
                            if( ptclObj_k.tagsDict["fftype"]  == 'O' ):
                                if( update_chr ): ptclObj_k.charge = q_O
                                if(set_resnames): ptclObj_k.tagsDict["resname"] = residue
                            if( ptclObj_k.tagsDict["fftype"]  == 'CT' ):
                                if( update_chr ): ptclObj_k.charge = q_BC
                                if(set_resnames): ptclObj_k.tagsDict["resname"] = residue
                                Nk_o = cov_nbindx[pid_k]
                                Nk_f = cov_nbindx[pid_k+1] - 1
                                for indx_k in range( Nk_o,Nk_f+1):
                                    pid_l = cov_nblist[indx_k]
                                    ptclObj_l = strucC.ptclC[pid_l]
                                    ELCNT_l = top.calc_elcnt(pid_l,strucC,cov_nblist,cov_nbindx)
                                    if( ptclObj_l.tagsDict["fftype"] == 'CT' ):
                                        #if( ELCNT[1] == 3 ):
                                        #    if( update_chr ): CHARGES[k] = q_CH3
                                        #    RESID[l] = residue
                                        if( ELCNT_l[1] == 2 ):
                                            if( update_chr ): ptclObj_l.charge = q_CH2
                                            if(set_resnames): ptclObj_l.tagsDict["resname"] = residue
                                            Nl_o = cov_nbindx[pid_l]
                                            Nl_f = cov_nbindx[pid_l+1] - 1
                                            for indx_m in range( Nl_o,Nl_f+1):
                                                pid_m = cov_nblist[indx_m]
                                                ptclObj_m = strucC.ptclC[pid_m]
                                                if( ptclObj_m.tagsDict["number"]  == 1 ):
                                                    #ATYPE[k] = 'H1'
                                                    #if( update_chr ): CHARGES[k] = q_OHMe
                                                    if(set_resnames):  ptclObj_m.tagsDict["resname"]  = residue
                                                    if( debug ):
                                                        print " changing " ,  ptclObj_l.tagsDict["number"],ptclObj_l.tagsDict["symbol"],ptclObj_l.tagsDict["fftype"],ptclObj_l.tagsDict["resname"]
                                                        print "    to ",residue
                                                        #sys.exit('set_pmmacharges')
                                    # Hack to not reasign terminal methyls 
                                    if(ptclObj_l.tagsDict["fftype"] == 'CT' and ptclObj_l.tagsDict["gtype"].strip() == methyl_C  ):
                                        #if( ELCNT[1] == 3 ):
                                        #    if( update_chr ): CHARGES[k] = q_CH3
                                        if( ELCNT_l[1] == 3 ):
                                            if( update_chr ): ptclObj_l.charge  = q_CH3
                                            if( debug ):
                                                print " found pmma methyl ",  ptclObj_l.tagsDict["number"],ptclObj_l.tagsDict["symbol"],ptclObj_l.tagsDict["fftype"],ptclObj_l.tagsDict["resname"]
                                            if(set_resnames): ptclObj_l.tagsDict["resname"] = residue
                                            Nl_o = cov_nbindx[pid_l]
                                            Nl_f = cov_nbindx[pid_l+1] - 1
                                            for indx_m in range( Nl_o,Nl_f+1):
                                                pid_m = cov_nblist[indx_m]
                                                ptclObj_m = strucC.ptclC[pid_m]
                                                if( ptclObj_m.tagsDict["number"]  == 1 ):
                                                    #ATYPE[k] = 'H1'
                                                    if( update_chr ): ptclObj_m.charge = q_H3
                                                    if(set_resnames): ptclObj_m.tagsDict["resname"]  = residue
                                                    if( debug ):
                                                        print " changing " , ptclObj_l.tagsDict["number"],ptclObj_l.tagsDict["symbol"],ptclObj_l.tagsDict["fftype"],ptclObj_l.tagsDict["resname"]
                                                        print "    to ",residue
                                                        #sys.exit('set_pmmacharges')



    # Hack
    for pid_i, ptclObj_i  in strucC.ptclC:
        if ( ptclObj_i.tagsDict["fftype"].strip() == methyl_C ):
            N_o = cov_nbindx[pid_i]
            N_f = cov_nbindx[pid_i+1] - 1
            for indx in range( N_o,N_f+1):
                pid_j = cov_nblist[indx]
                ptclObj_j = strucC.ptclC[pid_j]
                if (ptclObj_j.tagsDict["fftype"].strip()== 'C2' ):
                    if( update_chr ): ptclObj_i.charge = q_CH3
                    if(set_resnames):  ptclObj_j.tagsDict["resname"] = residue

    return strucC
                    

def biaryl_types(update_chr,strucC , ring_nblist, ring_nbindex,cov_nblist, cov_nbindx ):
    
    debug = 0

    #
    # label carbons 
    #
    
    if( debug): print "biaryl_types"
    n_rings = 0
    for pid_i, ptclObj_i  in strucC.ptclC:
        rnumb_i = ptclObj_i.tagsDict["ring"]
        if ( rnumb_i > n_rings ): n_rings = rnumb_i
    
    for r_numb in range( 1,n_rings+1 ):
        nring = top.calc_nnab(r_numb,ring_nbindex)
        ELCNT = top.calc_elcnt(r_numb,strucC,ring_nblist,ring_nbindex)
        
        if( debug ): print 'ring ',r_numb

        Nr_o = ring_nbindex[r_numb]
        Nr_f = ring_nbindex[r_numb+1] - 1
        for r_indx in range(Nr_o,Nr_f+1):
            pid_i = ring_nblist[r_indx]
            ptclObj_i =  strucC.ptclC[pid_i]
            NNAB_i = top.calc_nnab(pid_i,cov_nbindx)
            ELCNT_i = top.calc_elcnt(pid_i,strucC,cov_nblist,cov_nbindx)
            N_o = cov_nbindx[ pid_i ]
            N_f = cov_nbindx[ pid_i + 1 ] - 1
            NNAB_i_intra = 0
            for indx_j in range( N_o,N_f+1):
                pid_j = cov_nblist[indx_j]
                ptclObj_j =  strucC.ptclC[pid_j]
                if( ptclObj_j.tagsDict["ring"] == r_numb and ptclObj_j.tagsDict["number"] != 1 ):
                    NNAB_i_intra = NNAB_i_intra +  1
                    
            if( debug ): print " atom ",pid_i,ptclObj_i.tagsDict["fftype"],r_numb,NNAB_i,NNAB_i_intra
            if(  NNAB_i_intra == 2 ):
                # edge
                if( ptclObj_i.tagsDict["number"] == 6  and NNAB_i == 3 ):
                    if( ELCNT_i[16] == 1 and ELCNT_i[6] == 1 and ELCNT_i[1] == 1  ):
                        ptclObj_i.tagsDict["fftype"] = 'CP'
                    elif( ELCNT_i[7] == 1 and ELCNT_i[6] == 1 and ELCNT_i[1] == 1  ):
                        ptclObj_i.tagsDict["fftype"] = 'CW'
                    elif( ELCNT_i[7] == 1 and ELCNT_i[6] == 1  ):
                        ptclObj_i.tagsDict["fftype"] = 'CW'
                    elif( ELCNT_i[8] == 1 and ELCNT_i[6] == 1 and ELCNT_i[1] == 1  ):
                        ptclObj_i.tagsDict["fftype"] = 'CW'
                    elif( ELCNT_i[16] == 1 and ELCNT_i[6] == 2 ):
                        ptclObj_i.tagsDict["fftype"] = 'CP'
                    elif( ELCNT_i[6] == 2 ):
                        ptclObj_i.tagsDict["fftype"] = 'CA'
                if( ptclObj_i.tagsDict["number"] == 16  and NNAB_i == 2 and  ELCNT_i[6] == 2 ):
                    ptclObj_i.tagsDict["fftype"] = 'S'
                if( ptclObj_i.tagsDict["number"] == 7 and  NNAB_i == 3 and  ELCNT_i[6] == 3 ):
                    ptclObj_i.tagsDict["fftype"] = 'NS'
                if( ptclObj_i.tagsDict["number"] == 7 and  NNAB_i == 3 and  ELCNT_i[6] == 2 and  ELCNT_i[1] == 1 ):
                    ptclObj_i.tagsDict["fftype"] = 'NA'
                        
            if(  NNAB_i_intra == 3 ):
                # fussed 
                if( ptclObj_i.tagsDict["number"] == 6  and NNAB_i == 3 ):
                        ptclObj_i.tagsDict["fftype"] = 'CB'
                        
            if( debug ):
                print " atom ",atom_i,ptclObj_i.tagsDict["fftype"],NNAB_i,NNAB_i_intra
                        
        for r_indx in range(Nr_o,Nr_f+1):
            # relable secondary atoms 
            pid_i = ring_nblist[r_indx]
            ptclObj_i =  strucC.ptclC[pid_i]
            NNAB_i = top.calc_nnab(pid_i,cov_nbindx)
            ELCNT_i = top.calc_elcnt(pid_i,strucC,cov_nblist,cov_nbindx)
            N_o = cov_nbindx[ pid_i ]
            N_f = cov_nbindx[ pid_i + 1 ] - 1
            NNAB_i_intra = 0
            CP_cnt = 0 
            CW_cnt = 0 
            for indx_j in range( N_o,N_f+1):
                pid_j = cov_nblist[indx_j]
                ptclObj_j =  strucC.ptclC[pid_j]
                if( ptclObj_j.tagsDict["ring"] == r_numb and  ptclObj_j.tagsDict["number"] != 1  ):
                    NNAB_i_intra = NNAB_i_intra +  1
                    if( ptclObj_j.tagsDict["fftype"]== "CP" ): CP_cnt = CP_cnt  + 1
                    if( ptclObj_j.tagsDict["fftype"]== "CW" ): CW_cnt = CW_cnt  + 1
                    
            if(  NNAB_i_intra == 2 ):
                if( ptclObj_i.tagsDict["number"] == 6  and NNAB_i == 3 ):
                    if( CP_cnt > 0  or  CW_cnt > 0 ):
                        ptclObj_i.tagsDict["fftype"] = 'CS'
                        
            if( debug ):
                print "    relable atom ",atom_i,ptclObj_i.tagsDict["fftype"],NNAB_i,NNAB_i_intra,CP_cnt,CW_cnt
        
    debug = 0         
    if( debug ):
        sys.exit('atom_types.biaryl_types')
    return ( strucC )

def interring_types(update_chr,strucC , ring_nblist, ring_nbindex,cov_nblist, cov_nbindx ):
    """
    Set conjugated inter-ring carbons to type C! 
    """
      
    # Find Ring linkers
    debug = 0
    for pid_i, ptclObj_i  in strucC.ptclC:
        NNAB = top.calc_nnab(pid_i,cov_nbindx)
        ELCNT = top.calc_elcnt(pid_i,strucC,cov_nblist,cov_nbindx)
        N_o = cov_nbindx[ pid_i ]
        N_f = cov_nbindx[ pid_i + 1 ] - 1
        if ( ptclObj_i.tagsDict["number"] == 6 ) :
            if( NNAB == 3 ):  # if sp2
                for indx_j in range( N_o,N_f+1):
                    pid_j = cov_nblist[indx_j]
                    ptclObj_j = strucC.ptclC[pid_j]
                    if ( ptclObj_j.tagsDict["number"] == 6 )  :
                        NNAB_j = top.calc_nnab(pid_j,cov_nbindx)
                        if( NNAB_j == 3 ):  # if sp2
                            if(  ptclObj_i.tagsDict["ring"] !=  ptclObj_j.tagsDict["ring"] and ptclObj_i.tagsDict["ring"] != 0 and ptclObj_j.tagsDict["ring"] != 0): # ring linker
                                if(debug):
                                    print atom_i,' and ',pid_j,' link ',ptclObj_i.tagsDict["fftype"],ATYPE[pid_j]
                                    print ' ', ptclObj_i.tagsDict["ring"] , ptclObj_j.tagsDict["ring"]
                                ptclObj_i.tagsDict["fftype"] = 'C!'
                                ptclObj_j.tagsDict["fftype"]= 'C!'
                                
    return (strucC ) 

