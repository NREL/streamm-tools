#! /usr/bin/env python
"""
Asign atom types for force-field topologies 
"""

# Dr. Travis Kemper
# NREL
# Initial Date 05/07/2013
# travis.kemper@nrel.gov

def initialize_atype( ELN ):
    """
    Set atom type to atomic symbol for initialization purposes 
    """
    import elements

    debug = 0
    
    ATYPE= []
    elsymbol = elements.set_elsymbol()
    for atom_i in range( len(ELN) ):
        el_i = ELN[atom_i]
        atomic_symbol = elsymbol[el_i]
        ATYPE.append( atomic_symbol )
        if( debug): print atom_i,ATYPE[atom_i]
            
    return ATYPE

def oplsaa( update_chr,ELN,CHARGES,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB ):
    """
    Set OPLSaa atom types 
    """
    import sys,top
    
    ATYPE = []

    NA = len(ELN)

    ATYPE = initialize_atype( ELN )
    
    for atom_i in range(NA):
        NNAB = top.calc_nnab(atom_i,NBLIST,NBINDEX)
        ELCNT = top.calc_elcnt(atom_i,ELN,NBLIST,NBINDEX)
        if ELN[atom_i] == 6 :
           
            # simple guess based on coordination 
            if int(NNAB) == 4 :
                ATYPE[atom_i] = 'CT' # Alkane
                # refine guess based on nieghbors 
                if( ELCNT[1] == 4 ):                              # Methane 
                    ATYPE[atom_i] = 'CT'
                    if( update_chr ): CHARGES[atom_i] =  -0.24
                elif( ELCNT[6] == 1 and ELCNT[1] == 3):         # Methyl 
                    ATYPE[atom_i] = 'CT'
                    if( update_chr ): CHARGES[atom_i] =  -0.18
                elif(   ELCNT[6] == 2 and ELCNT[1] == 2 ):      
                    ATYPE[atom_i] = 'CT'
                    if( update_chr ): CHARGES[atom_i] =  -0.12
                elif( ELCNT[6] == 3 and ELCNT[1] == 1 ):         # Alkane 
                    ATYPE[atom_i] = 'CT'
                    if( update_chr ): CHARGES[atom_i] =  -0.06
                elif(  ELCNT[6] == 4 ):
                    ATYPE[atom_i] = 'CT'
                    if( update_chr ): CHARGES[atom_i] =  0.0
                elif(  ELCNT[1] == 1 and ELCNT[8] == 2 ):       #acetal
                    ATYPE[atom_i] = 'CO'
                    if( update_chr ): CHARGES[atom_i] =  -0.4

            if  int(NNAB) == 3 :
                r_numb = RING_NUMB[atom_i]
                if( r_numb != 0 ):
                    nring = top.ring_natoms(r_numb,RINGINDEX)
                    nring_mod = nring % 2     # to tell 6 adn 5 member rings apart
                    ATYPE[atom_i] = 'CA' # Aromatic C
                    if( ELCNT[6] == 2 and ELCNT[1] == 1 ):           # "Aromatic C"  
                        ATYPE[atom_i] = 'CA'
                        if( update_chr ): CHARGES[atom_i] =  -0.1150
                    if( ELCNT[6] == 3  and nring_mod == 0 ):            # "Naphthalene Fusion C"
                        ATYPE[atom_i] = 'CA'
                        if( update_chr ): CHARGES[atom_i] =  0.0
                    if( nring == 5 or nring == 8 or nring == 12 ):   # Ring with 5 membered parts
                        if( ELCNT[6] == 3 ):            # Fused 5-6 rings 
                            ATYPE[atom_i] = 'CB'
                            if( update_chr ): CHARGES[atom_i] =  0.0
                        if( ELCNT[7] == 1 or ELCNT[16] == 1  ):      # Pyrrole or thiophene 
                            if( ELCNT[6] == 1 and ELCNT[1] == 1 ):             #imidazole C2
                                ATYPE[atom_i] = 'CW'   
                                if( update_chr ): CHARGES[atom_i] =  0.22
                            if( ELCNT[6] == 1 and ELCNT[8] == 1 ):             # Not sure about this one
                                ATYPE[atom_i] = 'CW'   
                                if( update_chr ): CHARGES[atom_i] =  0.22
                        if(  ELCNT[7] == 2 ):             #imidazole C2
                            ATYPE[atom_i] = 'CR'   
                            if( update_chr ): CHARGES[atom_i] =  0.22
                        if(  ELCNT[7] == 1 and ELCNT[16] == 1 ):             #imidazole C2
                            ATYPE[atom_i] = 'CR'
                            if( update_chr ): CHARGES[atom_i] =  0.22

                    if( ELCNT[6] == 2 and  ELCNT[8] == 1  ):            # CTD strangeness ... 
                        ATYPE[atom_i] = 'CW'
                        if( update_chr ): CHARGES[atom_i] =  0.0
                        
                    debug = 0
                    if( debug  ):
                        print " r_numb ",nring,r_numb, nring_mod, ELCNT[6], ELCNT[7], ATYPE[atom_i]
                        # sys.exit(' debug ')
                    
                else: # not aromatic
                    if( ELCNT[6] == 2 and ELCNT[1] == 1 ):          # diene 
                        ATYPE[atom_i] = 'C='
                        if( update_chr ): CHARGES[atom_i] =  0.0
                    elif( ELCNT[6] == 2 and ELCNT[8] == 1 ):          # Benzophenone
                        ATYPE[atom_i] = 'C'
                        if( update_chr ): CHARGES[atom_i] =  0.7
                    elif( ELCNT[6] == 1 and ELCNT[8] == 2 ):
                        ATYPE[atom_i] = 'C' # pmma
                        if( update_chr ): CHARGES[atom_i] =  0.7

            if int(NNAB) == 2 :
                ATYPE[atom_i] = 'C:'    # Allene
                if( ELCNT[6] == 1 and ELCNT[7] == 1 ):   # "Benzonitrile -CN"  
                    ATYPE[atom_i] = 'CZ'
                    
            if int(NNAB) == 1 :
                ATYPE[atom_i] = '' # Aromatic C
                print " WARNING!!! carbon index ",atom_i," bonded to single atom "
        #
        # label oxygens
        #
        if( ELN[atom_i] == 8 ):
            if int(NNAB) == 1 :
                ATYPE[atom_i] = 'O' # double bonded
                if( update_chr ): CHARGES[atom_i] =  -0.5
            if int(NNAB) == 2 :
                ATYPE[atom_i] = 'OS' # ether
                if( update_chr ): CHARGES[atom_i] =  -0.5
                if( ELCNT[1] == 1 ):
                    ATYPE[atom_i] = 'OH' # Alcohol
                    if( update_chr ): CHARGES[atom_i] =   -0.6830
                if( ELCNT[8] == 1 ):
                    ATYPE[atom_i] = 'O2' # Carboxylate
                    if( update_chr ): CHARGES[atom_i] =   -0.800
                if( ELCNT[16] == 1 ):
                    ATYPE[atom_i] = 'OY' # Sulfoxide
                    if( update_chr ): CHARGES[atom_i] =   -0.4200
            if int(NNAB) == 3 :
                if( ELCNT[7] == 1  ):
                    ATYPE[atom_i] = 'ON'
                    if( update_chr ): CHARGES[atom_i] =  -0.118
                        
        #
        # label nitrogens 
        #
        if ELN[atom_i] == 7 :
            if int(NNAB) == 3 :      # amide
                ATYPE[atom_i] = 'N' 
                if( ELCNT[1] == 3 ): # Ammonia NH3"   
                    ATYPE[atom_i] = 'NT'
                if( ELCNT[1] == 2 ): # -NH2
                    ATYPE[atom_i] = 'N2'
                if( ELCNT[1] == 1 ): # -NH2
                    ATYPE[atom_i] = 'N2'
                    
            r_numb = RING_NUMB[atom_i]
            if( r_numb != 0 ):
                RN_o = RINGINDEX[r_numb]
                RN_f = RINGINDEX[r_numb+1] - 1
                nring = RN_f - RN_o + 1
                nring_mod = nring % 2     # to tell 6 adn 5 member rings apart
                if ( nring == 5 or nring == 8 or nring == 12 ):   # Ring with 5 membered parts
                    if ( int(NNAB) == 2 ):
                        ATYPE[atom_i] = 'NC'      # "Imidazole N3"     
                    if  int(NNAB) == 3 :
                        ATYPE[atom_i] = 'NA'      # "Imidazole N1"      
                else:
                    if ( int(NNAB) == 2  ):
                        ATYPE[atom_i] = 'NB'
            if ( int(NNAB) == 1 and ELCNT[6] == 1 ) :      # Nitrile 
                ATYPE[atom_i] = 'NZ'
                

        #
        # label sulfurs
        #
        if( ELN[atom_i] == 16 ):
            if int(NNAB) == 2 :
                ATYPE[atom_i] = 'S'   #  Thioether RSR (UA)
                if( update_chr ): CHARGES[atom_i] =-0.4350
                if( ELCNT[1] == 1  ):
                    ATYPE[atom_i] = 'SH'
                    if( update_chr ): CHARGES[atom_i] =  -0.335
                if( ELCNT[1] == 2  ):
                    ATYPE[atom_i] = 'SH'
                    if( update_chr ): CHARGES[atom_i] =  -0.470
                        
            if ( int(NNAB) == 4 and ELCNT[8] == 2 and ELCNT[7] == 1 ):
                ATYPE[atom_i] = 'SY'   #  "Sulfonamide -SO2N<" 
            if ( int(NNAB) == 3 and ELCNT[8] == 1 ):
                ATYPE[atom_i] = 'SZ'   #  Sulfoxide   
                           
        # Label chlorine 
        if ( ELN[atom_i] == 17):
            ATYPE[atom_i] = 'Cl'
            if( update_chr ): CHARGES[atom_i] =  -0.2

        # Label Fluorine   
        if ( ELN[atom_i] == 9):
            ATYPE[atom_i] = 'F'
            if( update_chr ): CHARGES[atom_i] =  -0.2057

    #
    # label hydrogens
    #
    for atom_i in range(NA):
        NNAB = top.calc_nnab(atom_i,NBLIST,NBINDEX)
        ELCNT = top.calc_elcnt(atom_i,ELN,NBLIST,NBINDEX)
        if( ELN[atom_i] == 1 ):
            if ( NNAB > 1  ):
                sys.exit(' over coordinated H')
            atom_j = NBLIST[ NBINDEX[atom_i] ]
            ELCNT_j =  top.calc_elcnt(atom_j,ELN,NBLIST,NBINDEX)
            if ( ATYPE[atom_j] == 'CA' ):
                ATYPE[atom_i] = 'HA' #
                if( update_chr ): CHARGES[atom_i] =  0.115
            if ( ATYPE[atom_j] == 'CT' ):
                ATYPE[atom_i] = 'HC' #
                if( update_chr ): CHARGES[atom_i] =  0.06
            if ( ATYPE[atom_j] == 'CW' ):
                ATYPE[atom_i] = 'HA' #
            if ( ATYPE[atom_j] == 'CS' ):
                ATYPE[atom_i] = 'HA' #
            if( ELCNT_j[6] == 2 and ELCNT_j[8] == 1  ): # "Ester -OCH<"
                ATYPE[atom_i] = 'H1'
                if( update_chr ): CHARGES[atom_j] =  0.03
            if ( ATYPE[atom_j] == 'SH' ):
                ATYPE[atom_i] = 'HS' #
                if( update_chr ): CHARGES[atom_i] = 0.1550
            if ( ELN[atom_j] == 7 ):
                ATYPE[atom_i] = 'H' #
            if ( ATYPE[atom_j] == 'NA' ):
                ATYPE[atom_i] = 'H' #
            if ( ATYPE[atom_j] == 'N2' ):
                ATYPE[atom_i] = 'H' #
            if ( ATYPE[atom_j] == 'OH' ):
                ATYPE[atom_i] = 'HO' #
                    
                    
    # relabel based on neighbors
    for atom_i in range(NA):
        NNAB = top.calc_nnab(atom_i,NBLIST,NBINDEX)
        ELCNT = top.calc_elcnt(atom_i,ELN,NBLIST,NBINDEX)
        N_o = NBINDEX[ atom_i ]
        N_f = NBINDEX[ atom_i + 1 ] - 1
        nring_i = RING_NUMB[atom_i]
        if ( ATYPE[atom_i] == 'CA' ):
            for indx_j in range( N_o,N_f+1):
                atom_j = NBLIST[indx_j]
                ELCNT_j =  top.calc_elcnt(atom_j,ELN,NBLIST,NBINDEX)
                nring_j = RING_NUMB[atom_j]
                if ( ATYPE[atom_j] == 'CW' ):
                    if( ELCNT[1] == 1 and ELCNT[6] == 2 ):
                        ATYPE[atom_i] = 'CS'
                if(  ELCNT[6] == 2 and nring_i == nring_j and ELN[atom_j] != 6 ):          # fussed 
                    if( debug):
                        print ' CB ',atom_i+1, ELCNT[6] ,  nring_i , nring_j 
                    ATYPE[atom_i] = 'CB'
            if ( ATYPE[atom_i] == 'CS' ):
                for indx_j in range( N_o,N_f+1):
                    atom_j = NBLIST[indx_j]
                    ELCNT_j = top.calc_elcnt(atom_j,ELN,NBLIST,NBINDEX)
                    if (  ATYPE[atom_j] == 'CA' ):
                        if ( ELCNT_j[6] == 3 ):
                            ATYPE[atom_j] = 'CS'
                        if ( ELCNT_j[6] == 2 and ELCNT_j[1] == 1 ):
                            ATYPE[atom_j] = 'CS'

        if ATYPE[atom_i] == 'CT' :
            for indx_j in range( N_o,N_f+1):
                atom_j = NBLIST[indx_j]
                if ( ATYPE[atom_j] == 'O' ):
                    if( ELCNT[1] == 3 and ELCNT[8] == 1 ):
                        ATYPE[atom_i] = 'C3' # Methyloxide
                        
                        
    # Find Ring linkers
    debug = 0
    for atom_i in range(NA):
        NNAB = top.calc_nnab(atom_i,NBLIST,NBINDEX)
        ELCNT = top.calc_elcnt(atom_i,ELN,NBLIST,NBINDEX)
        N_o = NBINDEX[ atom_i ]
        N_f = NBINDEX[ atom_i + 1 ] - 1
        if ( ATYPE[atom_i] == 'CA' or  ATYPE[atom_i] == 'CB'  or  ATYPE[atom_i] == 'CW' ) :
            if( NNAB == 3 ):  # if sp2
                for indx_j in range( N_o,N_f+1):
                    atom_j = NBLIST[indx_j]
                    if ( ATYPE[atom_j] == 'CA' or  ATYPE[atom_j] == 'CB'  or  ATYPE[atom_j] == 'CW' )  :
                        NNAB_j = top.calc_nnab(atom_j,NBLIST,NBINDEX)
                        if( NNAB_j == 3 ):  # if sp2
                            if(  RING_NUMB[atom_i] !=  RING_NUMB[atom_j] ): # ring linker
                                if(debug):
                                    print atom_i,' and ',atom_j,' link ',ATYPE[atom_i],ATYPE[atom_j]
                                    print ' ', RING_NUMB[atom_i] , RING_NUMB[atom_j]
                                ATYPE[atom_i] = 'C!'
                                ATYPE[atom_j] = 'C!'                                

    # relabel based on neighbors
    debug = 0
    for atom_i in range(NA):
        if ( ATYPE[atom_i] == 'C!' ):
            r_numb = RING_NUMB[atom_i]
            nring =  top.ring_natoms(r_numb,RINGINDEX)
            ELCNT = top.calc_elcnt(atom_i,ELN,NBLIST,NBINDEX)
            if( debug ):
                print ' linker found ',atom_i
                print ' nring,ELCNT[6],ELCNT ',nring,ELCNT[6],ELCNT 
            if( nring == 5 or nring == 8 or nring == 12 ):   # Ring with 5 membered parts
                if( ELCNT[7] == 1 or ELCNT[16] == 1  ):      # Pyrrole or thiophene 
                    if( ELCNT[6] == 2  ):             #imidazole C2
                        N_o = NBINDEX[ atom_i ]
                        N_f = NBINDEX[ atom_i + 1 ] - 1
                        for indx_j in range( N_o,N_f+1):
                            atom_j = NBLIST[indx_j]
                            if( ATYPE[atom_j] == 'CA' ):
                                if( debug ): print " changing ",ATYPE[atom_j] , " to CS"
                                ATYPE[atom_j] = 'CS'   
                                if( update_chr ): CHARGES[atom_i] =  0.22

    if(debug):   sys.exit('linkers')
  
    debug = 0
    if(debug):
        for atom_i in range(NA):
            print atom_i,ATYPE[atom_i],ELN[atom_i]
        sys.exit('debug')

    # Check for unidentified atoms
    for atom_i in range(NA):
        if ( ATYPE[atom_i] == '?' ):
            print ' atom ', atom_i , ELN[atom_i],' unknow '
            sys.exit(' Unknow atom ')
         
    return (ATYPE,CHARGES)

def set_ptmatypes( update_chr, ELN,ASYMB, ATYPE,GTYPE,RESID,CHARGES,AMASS,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB,CG_SET,CHARN ):
    """
    Set nitroxyl atom types for TEMPO 
    """
    
    import sys,top 
    import sys,top 
        
        
    residue = 'PTMA'

    debug = 0
    
    
    NA = len(ELN)

    # Set default charges, which will be over writen in some cases 
    for atom_i in range(NA):

        if(debug): print ASYMB[atom_i]
        if ATYPE[atom_i] == 'HC' :
            if( update_chr ): CHARGES[atom_i] =  0.06
        if ATYPE[atom_i] == 'CT' :
            if( update_chr ): CHARGES[atom_i] =  -0.18
        if ASYMB[atom_i].strip() == 'LP' :
            ATYPE[atom_i] = 'LP'
            RESID[atom_i] = residue
            AMASS[atom_i] = 9.0 
            if( update_chr ): CHARGES[atom_i] =  -0.11
        
    for atom_i in range(NA):
        NNAB = top.calc_nnab(atom_i,NBLIST,NBINDEX)
        ELCNT = top.calc_elcnt(atom_i,ELN,NBLIST,NBINDEX)
        #
        # label nitrogens 
        #
        if ( ELN[atom_i] == 7 ):
            if( ELCNT[8] == 1 and ELCNT[6] == 2  ):
                Ni_o = NBINDEX[ atom_i ]
                Ni_f = NBINDEX[ atom_i+1  ] - 1
                nitroxide = 0
                cation = 0
                for indx_i in range( Ni_o,Ni_f+1):
                    atom_j = NBLIST[indx_i]
                    if ( ELN[atom_j] == 8 ):
                        NNAB_j = top.calc_nnab(atom_j,NBLIST,NBINDEX)
                        ELCNT_j = top.calc_elcnt(atom_j,ELN,NBLIST,NBINDEX)
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

                    ATYPE[atom_i] = 'NN'
                    if( update_chr ): CHARGES[atom_i] = q_NN
                    RESID[atom_i] = residue
                    for indx_i in range( Ni_o,Ni_f+1):
                        atom_j = NBLIST[indx_i]
                        
                        if ( ELN[atom_j] == 8 ):
                            ATYPE[atom_j] = 'ON'
                            if( update_chr ): CHARGES[atom_j] = q_ON
                                    
                            RESID[atom_j] = residue
                        if ( ELN[atom_j] == 6 ):
                            ATYPE[atom_j] = 'CT'
                            if( update_chr ): CHARGES[atom_j] = q_CNN
                            RESID[atom_j] = residue
                            Nj_o = NBINDEX[ atom_j ]
                            Nj_f = NBINDEX[ atom_j +1  ] - 1
                            for indx_j in range(Nj_o,Nj_f+1):
                                atom_k = NBLIST[indx_j]
                                RESID[atom_k] = residue
                                ELCNT_k = top.calc_elcnt(atom_k,ELN,NBLIST,NBINDEX)
                                if(  ELN[atom_k] == 6 ):
                                    if ( ELCNT_k[1] == 3 ):
                                        if( update_chr ): CHARGES[atom_k] =  q_CTN
                                        RESID[atom_k] = residue
                                        Nk_o = NBINDEX[ atom_k ]
                                        Nk_f = NBINDEX[ atom_k+1  ] - 1
                                        for indx_k in range(Nk_o,Nk_f+1):
                                            atom_l = NBLIST[indx_k]
                                            if( ELN[atom_l] == 1):
                                                ATYPE[atom_l] = 'HN'
                                                if( update_chr ): CHARGES[atom_l] =  q_HN
                                                RESID[atom_l] = residue
                                            
                                    if ( ELCNT_k[1] == 2  and ELCNT_k[1] == 2 ):       
                                        if( update_chr ): CHARGES[atom_k] =  q_CT2
                                        RESID[atom_k] = residue
                                        Nk_o = NBINDEX[ atom_k ]
                                        Nk_f = NBINDEX[ atom_k+1  ] - 1
                                        for indx_k in range(Nk_o,Nk_f+1):
                                            atom_l = NBLIST[indx_k]
                                            RESID[atom_l] = residue
                                            if( ELN[atom_l] == 1 ):
                                                if( update_chr ): CHARGES[atom_l] =  q_HCT2
                                                RESID[atom_l] = residue
                                            if( ELN[atom_l] == 6 and atom_l != atom_k and atom_l != atom_j ):
                                                ATYPE[atom_l] = 'CT'
                                                if( update_chr ): CHARGES[atom_l] =  q_Cp
                                                RESID[atom_l] = residue
                                                Nl_o = NBINDEX[ atom_l ]
                                                Nl_f = NBINDEX[ atom_l+1  ] - 1
                                                for indx_l in range(Nl_o,Nl_f+1):
                                                    atom_m = NBLIST[indx_l]
                                                    RESID[atom_m] = residue
                                                    if( ELN[atom_m] == 1 ):
                                                        ATYPE[atom_m] = 'H1'
                                                        if( update_chr ): CHARGES[atom_m] =  q_H1
                                                        RESID[atom_m] = residue
                                                    if( ELN[atom_m] == 8 ):
                                                        ATYPE[atom_m] = 'OS'
                                                        if( update_chr ): CHARGES[atom_m] =  q_OS
                                                        RESID[atom_m] = residue
                                                        Nm_o = NBINDEX[ atom_m ]
                                                        Nm_f = NBINDEX[ atom_m+1  ] - 1
                                                        for indx_m in range(Nm_o,Nm_f+1):
                                                            atom_n = NBLIST[indx_m]
                                                            RESID[atom_n] = residue
                                                            if( ELN[atom_n] == 6 and atom_n != atom_l ):
                                                                ATYPE[atom_n] = 'C'
                                                                if( update_chr ): CHARGES[atom_n] =  q_C
                                                                RESID[atom_n] = residue
                                                                Nn_o = NBINDEX[ atom_n ]
                                                                Nn_f = NBINDEX[ atom_n+1  ] - 1
                                                                for indx_n in range(Nn_o,Nn_f+1):
                                                                    atom_o = NBLIST[indx_n]
                                                                    if( ELN[atom_o] == 8 and atom_o != atom_m ):
                                                                        ATYPE[atom_o] = 'O'
                                                                        if( update_chr ): CHARGES[atom_o] =  q_O
                                                                        RESID[atom_o] = residue
                                                                    if( ELN[atom_o] == 6 ):
                                                                        ATYPE[atom_o] = 'CT'
                                                                        if( update_chr ): CHARGES[atom_o] =  q_Cb1
                                                                        RESID[atom_o] = 'BCK' #residue
                                                                        No_o = NBINDEX[ atom_o ]
                                                                        No_f = NBINDEX[ atom_o+1  ] - 1
                                                                        for indx_o in range(No_o,No_f+1):
                                                                            atom_p = NBLIST[indx_o]
                                                                            ELCNT_p = top.calc_elcnt(atom_p,ELN,NBLIST,NBINDEX)
                                                                            
                                                                            
                                                                            
                                                                            if( ELN[atom_p] == 6 and atom_p != atom_n ):
                                                                        
                                                                                if( ELCNT_p[1] == 2 ):
                                                                                    ATYPE[atom_p] = 'CT'
                                                                                    if( update_chr ): CHARGES[atom_p] =  q_CH2
                                                                                    RESID[atom_p] = 'BCK' #residue
                                                                                    Np_o = NBINDEX[ atom_p ]
                                                                                    Np_f = NBINDEX[ atom_p+1  ] - 1
                                                                                    for indx_p in range(Np_o,Np_f+1):
                                                                                        atom_q = NBLIST[indx_p]
                                                                                        if( ELN[atom_q] == 1 ):
                                                                                            #ATYPE[atom_q] = 'H1'
                                                                                            #if( update_chr ): CHARGES[atom_q] = q_H3
                                                                                            RESID[atom_q] = 'BCK' # residue
                                                                                        
                                                                                if( GTYPE[atom_p].strip()  == options.methyl_C.strip() ): #or  # ELCNT_l[1] == 3 ):
                                                                                    
                                                                                    if( update_chr ): CHARGES[atom_p] = q_CH3
                                                                                    RESID[atom_p] = residue
                                                                                    Np_o = NBINDEX[atom_p]
                                                                                    Np_f = NBINDEX[atom_p+1] - 1
                                                                                    for indx_p in range( Np_o,Np_f+1):
                                                                                        atom_q = NBLIST[indx_p]
                                                                                        if( ELN[atom_q] == 1 ):
                                                                                            #ATYPE[atom_q] = 'H1'
                                                                                            if( update_chr ): CHARGES[atom_q] = q_H3
                                                                                            RESID[atom_q] = residue
                                                                                            
                                                    

    debug=0
    if( debug ):
        print NA
        for atom_i in range(len(ELN)):            
            print atom_i+1 ,ELN[atom_i],GTYPE[atom_i],ATYPE[atom_i],CHARGES[atom_i],RESID[atom_i]
        sys.exit('ptma_types ')


    #
    ## Set charge groups 
    #n_groups = 0
    #    
    ## set heavy atom groups
    #for i in range(NA):
    #    q_g = CHARN[atom_i]
    #    if [ q_g > n_groups ]: n_groups =  q_g 
    #if( options.verbose ):
    #    print " "
    ## set heavy atom groups
    #for i in range(NA):
    #    if( CG_SET[i] ):
    #        NNAB = calc_nnab(i,NBLIST,NBINDEX)
    #        ELCNT = calc_elcnt(i,ELN,NBLIST,NBINDEX)
    #        N_o = NBINDEX[i]
    #        N_f = NBINDEX[i+1] - 1
    #        if( ATYPE[i] == 'NN' ):  # Methyl 
    #            n_groups = CHARN[i]
    #            CG_SET[i] = 0
    #            CHARN[i] = n_groups
    #            for j_indx in range( N_o,N_f+1):
    #                j = NBLIST[j_indx]
    #                #if ( ASYMB[j] == 'H' ):
    #                CHARN[j] = n_groups
    #                CG_SET[j] = 0
    #                    
                        
    return ( ATYPE,RESID,CHARGES,CG_SET,CHARN )

def set_ptma_imps(NA,NBLIST, NBINDEX,ELN,ASYMB,IMPS,IMPTYPE_F):
    """
    Set nitroxyl improper dihedrals 
    """
    
    import sys,top 
    import sys
    #
    debug = 0
    #
    # Generate impropers from neighbor list
    #
    for i in range(NA):
        N_o = NBINDEX[ i ]
        N_f = NBINDEX[ i + 1 ] - 1
        NNAB = N_f - N_o + 1
        if ( ELN[i] == 8 and NNAB == 3 ):
            IMP_ATOMS = []
            make_imp = 0
            for j_indx in range( N_o,N_f+1):
                j = NBLIST[j_indx]
                if ( ASYMB[j].strip() == 'LP' ):
                    make_imp = 1
                    IMP_ATOMS.append(j)
            for j_indx in range( N_o,N_f+1):
                j = NBLIST[j_indx]
                if ( ASYMB[j].strip() == 'N' ):
                    make_imp = 1
                    IMP_ATOMS.append(j)
            if( make_imp ):
                IMPS.append([IMP_ATOMS[0],IMP_ATOMS[1],i,IMP_ATOMS[2] ])
                IMPTYPE_F.append(2)

    if(debug):
        print len(IMPS),' impropers found '
        for imp_indx in range( len(IMPS)):
            print IMPS[imp_indx][0],IMPS[imp_indx][1],IMPS[imp_indx][2],IMPS[imp_indx][3],IMPTYPE_F[imp_indx]
        sys.exit('set_ptma_imps')
        

    return ( IMPS,IMPTYPE_F )
        

def set_iontypes(charge_total,update_chr, ELN, ATYPE,GTYPE,RESID,CHARGES,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB ):
    """
    Set atom types for ions 
    """
    import sys,top 
    
    NA = len(ELN)
    
    debug = 0


    
    # Find F
    for atom_i in range(NA):
        NNAB = top.calc_nnab(atom_i,NBLIST,NBINDEX)
        ELCNT = top.calc_elcnt(atom_i,ELN,NBLIST,NBINDEX)
        # Fluorine 
        if ( ELN[atom_i] == 9 and NNAB == 0 ):
             ATYPE[atom_i] = 'F'
             if( update_chr ): CHARGES[atom_i] = -1.0
             RESID[atom_i] = "FLU"
             
        if ( ELN[atom_i] == 5 and NNAB == 4 and ELCNT[9] == 4 ):
            ATYPE[atom_i] = 'B4'
            if( update_chr ): CHARGES[atom_i] =  1.192
            RESID[atom_i] = "BF4"            
            N_o = NBINDEX[atom_i]
            N_f = NBINDEX[atom_i+1] - 1
            for indx_j in range( N_o,N_f+1):
                atom_j = NBLIST[indx_j]
                ATYPE[atom_j] = 'F4'
                if( update_chr ): CHARGES[atom_j] = -0.548
                RESID[atom_j] = "BF4"            


        if ( ELN[atom_i] == 17  and NNAB == 4 and ELCNT[8] == 4 ):
            ATYPE[atom_i] = 'CP'
            if( update_chr ): CHARGES[atom_i] = q_H3
            RESID[atom_i] = "ClO4"            
            N_o = NBINDEX[atom_i]
            N_f = NBINDEX[atom_i+1] - 1
            for indx_j in range( N_o,N_f+1):
                atom_j = NBLIST[indx_j]
                ATYPE[atom_j] = 'OP'
                if( update_chr ): CHARGES[atom_i] = q_H3
                RESID[atom_j] = "ClO4"            


    debug = 0
    if(debug):
        for atom_i in range(NA):
            print atom_i,ATYPE[atom_i],ELN[atom_i],RESID[atom_i],CHARGES[atom_i]
        sys.exit('debug')

            
    return ( ATYPE,RESID,CHARGES )

def set_pmmatypes(update_chr, ELN, ATYPE,GTYPE,RESID,CHARGES,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB ):
    """
    Set atom types for PMMA monomer
    """
    import sys,top 
    
    NA = len(ELN)
    
    debug = 0
    update_chr = 1 

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

    # Find pmma segments and label charges
    for i in range(NA):
        atom_i = i 
        NNAB = top.calc_nnab(atom_i,NBLIST,NBINDEX)
        ELCNT = top.calc_elcnt(atom_i,ELN,NBLIST,NBINDEX)
        N_o = NBINDEX[i]
        N_f = NBINDEX[i+1] - 1
        if ATYPE[i] == 'OS' :
            pmma_oxygen = 0
            for indx in range( N_o,N_f+1):
                j = NBLIST[indx]
                atom_j = j
                NNAB_j =  top.calc_nnab(atom_j,NBLIST,NBINDEX)
                ELCNT_j = top.calc_elcnt(atom_j,ELN,NBLIST,NBINDEX)
                if( ELCNT_j[1] == 3 ):
                    pmma_oxygen = 1
            if( pmma_oxygen == 1 ):
                if( update_chr ): CHARGES[i] = q_OS
                RESID[i] = residue
                for indx in range( N_o,N_f+1):
                    j = NBLIST[indx]
                    if( ATYPE[j] == 'CT' ):
                        if( update_chr ): CHARGES[j] = q_OMe
                        RESID[j] = residue
                        Nj_o = NBINDEX[j]
                        Nj_f = NBINDEX[j+1] - 1
                        for indx_j in range( Nj_o,Nj_f+1):
                            k = NBLIST[indx_j]
                            if( ATYPE[k] == 'HC' ):
                                ATYPE[k] = 'H1'
                                if( update_chr ): CHARGES[k] = q_OHMe
                                RESID[k] = residue
                                #
 
                    if( ATYPE[j] == 'C' ):
                        if( update_chr ): CHARGES[j] = q_C
                        RESID[j] = residue
                        Nj_o = NBINDEX[j]
                        Nj_f = NBINDEX[j+1] - 1
                        for indx_j in range( Nj_o,Nj_f+1):
                            k = NBLIST[indx_j]
                            if( ATYPE[k] == 'O' ):
                                if( update_chr ): CHARGES[k] = q_O
                                RESID[k] = residue
                            if( ATYPE[k] == 'CT' ):
                                if( update_chr ): CHARGES[k] = q_BC
                                RESID[k] = residue
                                Nk_o = NBINDEX[k]
                                Nk_f = NBINDEX[k+1] - 1
                                for indx_k in range( Nk_o,Nk_f+1):
                                    l = NBLIST[indx_k]
                                    atom_l = l
                                    ELCNT_l = top.calc_elcnt(atom_l,ELN,NBLIST,NBINDEX)
                                    if( ATYPE[l] == 'CT' ):
                                        #if( ELCNT[1] == 3 ):
                                        #    if( update_chr ): CHARGES[k] = q_CH3
                                        #    RESID[l] = residue
                                        if( ELCNT_l[1] == 2 ):
                                            if( update_chr ): CHARGES[l] = q_CH2
                                            RESID[l] = residue
                                            Nl_o = NBINDEX[l]
                                            Nl_f = NBINDEX[l+1] - 1
                                            for indx_m in range( Nl_o,Nl_f+1):
                                                m = NBLIST[indx_m]
                                                if( ELN[m] == 1 ):
                                                    #ATYPE[k] = 'H1'
                                                    #if( update_chr ): CHARGES[k] = q_OHMe
                                                    RESID[m] = residue
                                                    if( debug ):
                                                        print " changing " , ATYPE[m],ELN[m],RESID[m],GTYPE[m],k
                                                        print "    to ",residue
                                                        sys.exit('set_pmmacharges')
                                    # Hack to not reasign terminal methyls 
                                    if( ATYPE[l] == 'CT' and GTYPE[l].strip() == 'C19' ):
                                        #if( ELCNT[1] == 3 ):
                                        #    if( update_chr ): CHARGES[k] = q_CH3
                                        if( ELCNT_l[1] == 3 ):
                                            if( update_chr ): CHARGES[l] = q_CH3
                                            RESID[l] = residue
                                            Nl_o = NBINDEX[l]
                                            Nl_f = NBINDEX[l+1] - 1
                                            for indx_m in range( Nl_o,Nl_f+1):
                                                m = NBLIST[indx_m]
                                                if( ELN[m] == 1 ):
                                                    #ATYPE[k] = 'H1'
                                                    if( update_chr ): CHARGES[m] = q_H3
                                                    RESID[m] = residue
                                                    if( debug ):
                                                        print " changing " , ATYPE[m],ELN[m],RESID[m],GTYPE[m],k
                                                        print "    to ",residue
                                                        sys.exit('set_pmmacharges')



    # Hack
    for i in range(NA):
        if ( GTYPE[i].strip() == 'C19' ):
            N_o = NBINDEX[i]
            N_f = NBINDEX[i+1] - 1
            for indx in range( N_o,N_f+1):
                j = NBLIST[indx]
                if ( GTYPE[j].strip() == 'C4' ):
                    #if( update_chr ): CHARGES[i] = q_CH3
                    RESID[i] = residue

    return ( ATYPE,RESID,CHARGES)
                    

def biaryl_types(update_chr, ATYPE, ELN,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB, CHARGES ):
    import sys,top
    
    NA = len(ELN)
    
    debug = 0

    #
    # label carbons 
    #
    
    if( debug): print "biaryl_types"
    n_rings = 0
    for atom_i in range(NA):
        rnumb_i = RING_NUMB[atom_i]
        if ( rnumb_i > n_rings ): n_rings = rnumb_i
    
    for r_numb in range( 1,n_rings+1 ):
        nring = top.ring_natoms(r_numb,RINGINDEX)
        ELCNT = top.calc_elcnt(r_numb,ELN,RINGLIST,RINGINDEX)
        
        if( debug ): print 'ring ',r_numb

        Nr_o = RINGINDEX[r_numb]
        Nr_f = RINGINDEX[r_numb+1] - 1
        for r_indx in range(Nr_o,Nr_f+1):
            atom_i = RINGLIST[r_indx]
            NNAB_i = top.calc_nnab(atom_i,NBLIST,NBINDEX)
            ELCNT_i = top.calc_elcnt(atom_i,ELN,NBLIST,NBINDEX)
            N_o = NBINDEX[ atom_i ]
            N_f = NBINDEX[ atom_i + 1 ] - 1
            NNAB_i_intra = 0
            for indx_j in range( N_o,N_f+1):
                atom_j = NBLIST[indx_j]
                if( RING_NUMB[atom_j] == r_numb and ELN[atom_j] != 1 ):
                    NNAB_i_intra = NNAB_i_intra +  1
                    
            if( debug ): print " atom ",atom_i,ATYPE[atom_i],r_numb,NNAB_i,NNAB_i_intra
            if(  NNAB_i_intra == 2 ):
                # edge
                if( ELN[atom_i] == 6  and NNAB_i == 3 ):
                    if( ELCNT_i[16] == 1 and ELCNT_i[6] == 1 and ELCNT_i[1] == 1  ):
                        ATYPE[atom_i] = 'CP'
                    elif( ELCNT_i[7] == 1 and ELCNT_i[6] == 1 and ELCNT_i[1] == 1  ):
                        ATYPE[atom_i] = 'CW'
                    elif( ELCNT_i[7] == 1 and ELCNT_i[6] == 1  ):
                        ATYPE[atom_i] = 'CW'
                    elif( ELCNT_i[8] == 1 and ELCNT_i[6] == 1 and ELCNT_i[1] == 1  ):
                        ATYPE[atom_i] = 'CW'
                    elif( ELCNT_i[16] == 1 and ELCNT_i[6] == 2 ):
                        ATYPE[atom_i] = 'CP'
                    elif( ELCNT_i[6] == 2 ):
                        ATYPE[atom_i] = 'CA'
                if( ELN[atom_i] == 16  and NNAB_i == 2 and  ELCNT_i[6] == 2 ):
                    ATYPE[atom_i] = 'S'
                if( ELN[atom_i] == 7 and  NNAB_i == 3 and  ELCNT_i[6] == 3 ):
                    ATYPE[atom_i] = 'NS'
                if( ELN[atom_i] == 7 and  NNAB_i == 3 and  ELCNT_i[6] == 2 and  ELCNT_i[1] == 1 ):
                    ATYPE[atom_i] = 'NA'
                        
            if(  NNAB_i_intra == 3 ):
                # fussed 
                if( ELN[atom_i] == 6  and NNAB_i == 3 ):
                        ATYPE[atom_i] = 'CB'
                        
            if( debug ):
                print " atom ",atom_i,ATYPE[atom_i],NNAB_i,NNAB_i_intra
                        
        for r_indx in range(Nr_o,Nr_f+1):
            # relable secondary atoms 
            atom_i = RINGLIST[r_indx]
            NNAB_i = top.calc_nnab(atom_i,NBLIST,NBINDEX)
            ELCNT_i = top.calc_elcnt(atom_i,ELN,NBLIST,NBINDEX)
            N_o = NBINDEX[ atom_i ]
            N_f = NBINDEX[ atom_i + 1 ] - 1
            NNAB_i_intra = 0
            CP_cnt = 0 
            CW_cnt = 0 
            for indx_j in range( N_o,N_f+1):
                atom_j = NBLIST[indx_j]
                if( RING_NUMB[atom_j] == r_numb and  ELN[atom_j] != 1  ):
                    NNAB_i_intra = NNAB_i_intra +  1
                    if( ATYPE[atom_j] == "CP" ): CP_cnt = CP_cnt  + 1
                    if( ATYPE[atom_j] == "CW" ): CW_cnt = CW_cnt  + 1
                    
            if(  NNAB_i_intra == 2 ):
                if( ELN[atom_i] == 6  and NNAB_i == 3 ):
                    if( CP_cnt > 0  or  CW_cnt > 0 ):
                        ATYPE[atom_i] = 'CS'
                        
            if( debug ):
                print "    relable atom ",atom_i,ATYPE[atom_i],NNAB_i,NNAB_i_intra,CP_cnt,CW_cnt
        
    debug = 0         
    if( debug ):
        sys.exit('atom_types.biaryl_types')
    return ( ATYPE , CHARGES )

def interring_types(update_chr, ATYPE, ELN,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB, CHARGES ):
    """
    Set conjugated inter-ring carbons to type C! 
    """
    import sys,top
      
    NA = len(ELN)
    
    # Find Ring linkers
    debug = 0
    for atom_i in range(NA):
        NNAB = top.calc_nnab(atom_i,NBLIST,NBINDEX)
        ELCNT = top.calc_elcnt(atom_i,ELN,NBLIST,NBINDEX)
        N_o = NBINDEX[ atom_i ]
        N_f = NBINDEX[ atom_i + 1 ] - 1
        if ( ELN[atom_i] == 6 ) :
            if( NNAB == 3 ):  # if sp2
                for indx_j in range( N_o,N_f+1):
                    atom_j = NBLIST[indx_j]
                    if ( ELN[atom_i] == 6 )  :
                        NNAB_j = top.calc_nnab(atom_j,NBLIST,NBINDEX)
                        if( NNAB_j == 3 ):  # if sp2
                            if(  RING_NUMB[atom_i] !=  RING_NUMB[atom_j] and RING_NUMB[atom_i] != 0 and RING_NUMB[atom_j] != 0): # ring linker
                                if(debug):
                                    print atom_i,' and ',atom_j,' link ',ATYPE[atom_i],ATYPE[atom_j]
                                    print ' ', RING_NUMB[atom_i] , RING_NUMB[atom_j]
                                ATYPE[atom_i] = 'C!'
                                ATYPE[atom_j] = 'C!'
                                
    return ( ATYPE , CHARGES ) 

