import mpiBase
import sys 
from simulationGaussian import SimulationGaussian
from buildingblocks import Buildingblock


def get_options():
    import os, os.path
    from optparse import OptionParser
 
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    parser.add_option("--in_fchk", dest="in_fchk", type="string", default="", help="Input .fchk GAUSSIAN file")
    parser.add_option("--out_id", dest="out_id", type="string", default="out", help="id")
    
    (options, args) = parser.parse_args()
        
    return options, args


def calc_nnab(i,NBINDEX):
    """
    Return number of nieghbors for a given atom 
    """
    #
    # Find number of elements 
    #
    N_o = NBINDEX[ i  ]
    N_f = NBINDEX[  i+1  ] - 1 
    NNAB = N_f - N_o + 1
    return NNAB


def calc_elcnt(i,strucC,NBLIST,NBINDEX):
    """
    Return 
    """
    
    import numpy
    #
    # Find number of elements 
    #
    ELCNT = numpy.zeros(120, dtype =int )
    N_o = NBINDEX[ i  ]
    N_f = NBINDEX[  i+1  ] - 1 
    for indx in range( N_o,N_f+1):
        j = NBLIST[indx]
        el_j = int( strucC.ptclC[j].tagsDict["number"] )
	if( el_j >= 0 ):
	    ELCNT[el_j] = ELCNT[el_j] + 1

    return ELCNT 

def oplsaa_atomtypes(strucC,cov_nblist, cov_nbindx):
    """
    Guess OPLS-aa atom types based on coordination and nearest nieghbors 
    """
    
    for pid_i, ptclObj_i  in strucC.ptclC:
        NNAB = calc_nnab(pid_i,cov_nbindx)
        ELCNT = calc_elcnt(pid_i,strucC,cov_nblist,cov_nbindx)
        #
        # label carbons 
        #
        if ptclObj_i.tagsDict["number"] == 6 :
            if int(NNAB) == 4 :
                ptclObj_i.tagsDict["fftype"] = 'CT' # Alkane
            if  int(NNAB) == 3 :
                ptclObj_i.tagsDict["fftype"] = 'CA'  # Conjugated 
            if int(NNAB) == 2 :
                ptclObj_i.tagsDict["fftype"] = 'C:'   # Allene

                    
            if int(NNAB) == 1 :
                ptclObj_i.tagsDict["fftype"] = '' # Aromatic C
                error_line =  " WARNING!!! carbon index ",pid_i," bonded to single atom "
                sys.exit(error_line)
        #
        # label oxygens
        #
        if( ptclObj_i.tagsDict["number"] == 8 ):
            if int(NNAB) == 1 :
                ptclObj_i.tagsDict["fftype"] = 'O' # double bonded
            if int(NNAB) == 2 :
                ptclObj_i.tagsDict["fftype"] = 'OS' # ether

        #
        # label nitrogens 
        #
        if ptclObj_i.tagsDict["number"] == 7 :
            if int(NNAB) == 3 :      # amide
                ptclObj_i.tagsDict["fftype"] = 'N' 

        #
        # label sulfurs
        #
        if( ptclObj_i.tagsDict["number"] == 16 ):
            if int(NNAB) == 2 :
                ptclObj_i.tagsDict["fftype"] = 'S'   #  Thioether RSR (UA)


    #
    # label hydrogens
    #
    for pid_i, ptclObj_i  in strucC.ptclC:
        NNAB = calc_nnab(pid_i,cov_nbindx)
        ELCNT = calc_elcnt(pid_i,strucC,cov_nblist,cov_nbindx)
        if( ptclObj_i.tagsDict["number"] == 1 ):
            if ( NNAB > 1  ):
                sys.exit(' over coordinated H')
            pid_j = cov_nblist[ cov_nbindx[pid_i] ]
            ptclObj_j = strucC.ptclC[pid_j]
            ELCNT_j =  calc_elcnt(pid_j,strucC,cov_nblist,cov_nbindx)

            if ( ptclObj_j.tagsDict["fftype"]== 'CA' ):
                ptclObj_i.tagsDict["fftype"] = 'HA' #
            if ( ptclObj_j.tagsDict["fftype"]== 'CT' ):
                ptclObj_i.tagsDict["fftype"] = 'HC' #


    # Check for unidentified atoms
    for pid_i, ptclObj_i  in strucC.ptclC:
        if ( ptclObj_i.tagsDict["fftype"] == '?' ):
            print ' atom ', pid_i , ptclObj_i.tagsDict["number"],' unknow '
            sys.exit(' Unknow atom ')

    return


def id_ring3(pid_o,ptclObj_o,strucC,cov_nblist, cov_nbindx):
    """
    Find atoms in conjugated rings 
    """
    
    max_paths = 1000
    
    debug = False
    
    R_SET = []
    RING_ATOMS = []
    BAD_PATH = []

    #a_i = pid_o -1 
    one = 1
    zero = 0
    atoms_in_ring = 0
    # relabel based on neighbors
    n_ptcl = len(strucC.ptclC)

    if(debug):
        print " for a system of %d partiles checking particle %d "%(n_ptcl,pid_o)
        
    for pid_i in range(n_ptcl+1): 
        R_SET.append(one)
        
    R_SET[pid_o] = 0 
    r_term = 1
    cnt = 0
    p_cnt = 0
    if(debug): print ' initializing ring ',pid_o," w NN ", calc_nnab(pid_o,cov_nbindx)

    last_i = pid_o
    NNAB_last = calc_nnab(pid_i,cov_nbindx)
    while ( r_term ):
        N_o = cov_nbindx[last_i]
        N_f = cov_nbindx[last_i+1] - 1

        if( debug): print " checking neighbors ",N_o,N_f

        for n_indx in range( N_o,N_f+1):
            j = cov_nblist[n_indx]
            cnt = cnt + 1
            # If non of the neighboring atoms of atom j are conjugated
            #    make the path as bad 
            if( cnt > NNAB_last+1 ):
                p_cnt += 1
                if(debug): print ' bad path found resetting at atom ',j,strucC.ptclC[j].tagsDict["number"]
                for ring_a in range(len(RING_ATOMS)):
                    j = RING_ATOMS[ring_a]
                    
                BAD_PATH.append(j)
                # Reset lists and counts
                RING_ATOMS = []
                for pid_k in range(n_ptcl):
                    R_SET[pid_k] = 1
                atoms_in_ring = 0
                R_SET[pid_o] = 0 
                cnt = 0
                last_i = pid_o
                
                if(debug): print '  resetting last_i to ',pid_o,ptclObj_o.tagsDict["number"]
                    
                for bad_i in range( len(BAD_PATH)):
                    j = BAD_PATH[bad_i]
                    R_SET[j] = 0
                    if(debug): print '  bad path atoms ',strucC.ptclC[j].tagsDict["number"]
                break

            # If path traces back to original atom ring has been found 
            if( atoms_in_ring  > 1 and j == pid_o ):
                r_term = 0
                if(debug): print ' ring found with ',atoms_in_ring 
                break

            number_j = strucC.ptclC[j].tagsDict["number"]
            NNAB_j = calc_nnab(j,cov_nbindx)
            ELCNT_j = calc_elcnt(j,strucC,cov_nblist,cov_nbindx)
            ring_type = False 
	    
            if(debug):
		print '           with ', number_j
		print '           with ', NNAB_j,
		print '           with ', R_SET[j]

            # Test if atom j is conjugated 
            if ( number_j == 6 and NNAB_j == 3 and R_SET[j] == 1 ):                      ring_type = True
            if ( number_j == 16 and NNAB_j == 2 and R_SET[j] == 1 and ELCNT_j[6] == 2 ): ring_type = True
            if ( number_j == 7 and NNAB_j >= 2 and R_SET[j] == 1 ):                      ring_type = True
            if( ring_type ):
                atoms_in_ring = atoms_in_ring + 1
                RING_ATOMS.append(j)
                R_SET[j] = 0
                last_i = j
                cnt = 0
                NNAB_last = NNAB_j 
                if(debug): print '   atom ',j,number_j,' added '
                break
            
        if (p_cnt > max_paths ):
            if(debug): print '     max paths considered ',max_paths
            r_term = 0
            
    return RING_ATOMS

def ring_type(pid_i,ptclObj_o,strucC,cov_nblist, cov_nbindx):
    """
    Criteria for a conjugated particle 
    """
    debug = False
    
    ring_type = False 
    number_j = strucC.ptclC[pid_i].tagsDict["number"]
    NNAB_j = calc_nnab(pid_i,cov_nbindx)
    ELCNT_j = calc_elcnt(pid_i,strucC,cov_nblist,cov_nbindx)

    if(debug):
        print '            atomic # %d  NN  %d '%(number_j, NNAB_j)

    # Test if atom j is conjugated 
    if ( number_j == 6 and NNAB_j == 3 ):                      ring_type = True
    if ( number_j == 16 and NNAB_j == 2 and ELCNT_j[6] == 2 ): ring_type = True
    if ( number_j == 7 and NNAB_j >= 2  ):                      ring_type = True

    return ring_type

def find_rings(strucC,cov_nblist, cov_nbindx):
    """
    Find conjugate rings
    """
    import sys

    debug = False
    
    RINGLIST = []
    RINGINDEX = []
    RING_NUMB = []
    
    IN_RING = []

    n_ptcl = len(strucC.ptclC)
    one = 0 
    zero = 0 
    ring_cnt = 0
    a_cnt = 0
    # relabel based on neighbors
    for pid_i, ptclObj_i  in strucC.ptclC:
        IN_RING.append(one)
        RINGLIST.append(zero)
        RING_NUMB.append(zero)
        RINGINDEX.append(zero)

    RING_NUMB.append(zero)
    RINGLIST.append(zero)
    RINGINDEX.append(zero)

        
    # relabel based on neighbors
    for pid_i, ptclObj_i  in strucC.ptclC:
        RING_ATOMS = []
        # If particle i is conjugated 
        if( ring_type(pid_i,ptclObj_i,strucC,cov_nblist, cov_nbindx) ):
            RING_ATOMS = id_ring3(pid_i,ptclObj_i,strucC,cov_nblist, cov_nbindx)
        if( debug ):
            print pid_i,RING_ATOMS
            
        if ( len(RING_ATOMS) > 1 ):
            # If member of ring alread part of another ring add new ring to existing ring
            r_numb = 0 
            for ring_a in range(len(RING_ATOMS)):
                j = RING_ATOMS[ring_a]
                if( RING_NUMB[j] != 0 ):
                    r_numb = RING_NUMB[j]
            if( r_numb == 0 ):
                ring_cnt = ring_cnt + 1
                r_numb = ring_cnt 
                #print ' new ring ',ring_cnt
                 
            for ring_a in range(len(RING_ATOMS)):
                j = RING_ATOMS[ring_a]
                if( RING_NUMB[j] == 0 ): 
                    a_cnt = a_cnt + 1
                    RING_NUMB[j] = r_numb

    a_cnt = 0
    for r_numb in range(1,ring_cnt+1):
        RINGINDEX[r_numb] = a_cnt + 1
        for i in range(n_ptcl):
            if( RING_NUMB[i] == r_numb ):
                a_cnt = a_cnt + 1
                RINGLIST[a_cnt] = i
                
    RINGINDEX[ring_cnt+1] = a_cnt + 1


    # Set attached H to have same ring #
    attach_H = False
    if( attach_H ):

        for r_numb in range(1,ring_cnt+1):
            N_o = RINGINDEX[r_numb]
            N_f = RINGINDEX[r_numb+1] - 1
            print ' Ring ', r_numb
            for r_indx in range(N_o,N_f+1):
                i = RINGLIST[r_indx]
                N_o = cov_nbindx[i]
                N_f = cov_nbindx[i+1] - 1
                for n_indx in range( N_o,N_f+1):
                    j = cov_nblist[n_indx]
                    if( strucC.ptclC[j].tagsDict["number"] == 1 ):
                        RING_NUMB[j] = r_numb
		

    # update particles 
    for pid_i, ptclObj_i  in strucC.ptclC:
        ptclObj_i.tagsDict["ring"] =  RING_NUMB[pid_i]
        if( debug ):
            print ' atom  ',pid_i,ptclObj_i.tagsDict["number"],ptclObj_i.tagsDict["ring"]
	    		
	    

        #debug = True
    if(debug):
        print ring_cnt , " rings found "
        for r_numb in range(1,ring_cnt+1):
            N_o = RINGINDEX[r_numb]
            N_f = RINGINDEX[r_numb+1] - 1
            print ' Ring ', r_numb
            for r_indx in range(N_o,N_f+1):
                i = RINGLIST[r_indx]
                print ' type ',strucC.ptclC[i].tagsDict["number"],' or ',i,RING_NUMB[i]
                #sys.exit('find_rings')
                
        for pid_i, ptclObj_i  in strucC.ptclC:    
            print ' atom  ',pid_i,ptclObj_i.tagsDict["number"],ptclObj_i.tagsDict["ring"]
	    
	sys.exit('top.find_rings')

    return ( RINGLIST, RINGINDEX  )


def initialize_fftags(strucC):
    """
    Add in particle properties need for force field 
    """
    for pid, pt_i  in strucC.ptclC:
        add_dict = pt_i.tagsDict
        add_dict["residue"] = 1
        add_dict["resname"] = "RES"
        add_dict["fftype"] = "??"
        add_dict["lmptype"] = -1
        add_dict["qgroup"] = 1
        add_dict["ffmass"] = pt_i.tagsDict["mass"] 
        add_dict["chain"] = 1
        add_dict["gtype"] = pt_i.tagsDict["symbol"]
        add_dict["linkid"] = ""
        pt_i.setTagsDict(add_dict)                    

    return 



def set_chargegroups(strucC , ring_nblist, ring_nbindex,cov_nblist, cov_nbindx ):
    """
    Set charge groups 
    """
    debug = False
    verbose = True 
    
    if( verbose ):
	print "   Setting charge groups "
    
    
    # initialize 
    n_groups = 0


    #CG_SET = []
    #for pid_i, ptclObj_i  in strucC.ptclC:
    #        CG_SET.append(True) 

    for pid_i, ptclObj_i  in strucC.ptclC:
        ptclObj_i.tagsDict["qgroup"] = 0
        
    # set rings as charge groups
    for pid_i, ptclObj_i  in strucC.ptclC:
        # If unset set chrge group 
        #    if( CG_SET[i] ):
        if( ptclObj_i.tagsDict["qgroup"] == 0  ):
            if( ptclObj_i.tagsDict["ring"] > 0 ):
                ptclObj_i.tagsDict["qgroup"] = ptclObj_i.tagsDict["ring"]
                if( debug ):
                    print " setting qgroup  ring ",pid_i,ptclObj_i.tagsDict["symbol"],ptclObj_i.tagsDict["qgroup"],ptclObj_i.tagsDict["ring"]
                if(  ptclObj_i.tagsDict["ring"] >  n_groups ):
                    n_groups = n_groups + 1
                
    #if( debug): sys.exit('ring charge groups')
            
    # set methyls and alkyls 
    for pid_i, ptclObj_i  in strucC.ptclC:
        # If unset set chrge group
        #    if( CG_SET[i] ):
        if( ptclObj_i.tagsDict["qgroup"] == 0  ):
            NNAB = calc_nnab(pid_i,cov_nbindx)
            ELCNT = calc_elcnt(pid_i,strucC,cov_nblist,cov_nbindx)
            N_o = cov_nbindx[pid_i]
            N_f = cov_nbindx[pid_i+1] - 1
            if( ptclObj_i.tagsDict["number"] == 6 and ELCNT[1] > 0 ):  # CH_n ( Methyls / alkyl ) 
                n_groups = n_groups + 1
                #CG_SET[i] = 0
                ptclObj_i.tagsDict["qgroup"] = n_groups
                for j_indx in range( N_o,N_f+1):
                    pid_j = cov_nblist[j_indx]
                    ptclObj_j = strucC.ptclC[pid_j]
                    if (ptclObj_j.tagsDict["number"] == 1 ):
                        ptclObj_j.tagsDict["qgroup"] = n_groups
                        #CG_SET[j] = 0
			
            if(  ptclObj_i.tagsDict["fftype"]  == 'C' ):  # 
                n_groups = n_groups + 1
                #CG_SET[i] = 0
                ptclObj_i.tagsDict["qgroup"] = n_groups
                for j_indx in range( N_o,N_f+1):
                    pid_j = cov_nblist[j_indx]
                    ptclObj_j = strucC.ptclC[pid_j]
                    if (ptclObj_j.tagsDict["number"] == 8 ):
                        ptclObj_j.tagsDict["qgroup"] = n_groups

                    #if ( ELN[j] == 8 ):
                    #    CHARN[j] = n_groups
                    #    CG_SET[j] = 0

            if(  ptclObj_i.tagsDict["fftype"] == 'ON' ):  # nitroxide
                n_groups = n_groups + 1
                #CG_SET[i] = 0
                ptclObj_i.tagsDict["qgroup"] = n_groups
                for j_indx in range( N_o,N_f+1):
                    pid_j = cov_nblist[j_indx]
                    ptclObj_j = strucC.ptclC[pid_j]
                    if (ptclObj_j.tagsDict["qgroup"] == 0 ):
                        ptclObj_j.tagsDict["qgroup"] = n_groups

                    #j = cov_nblist[j_indx]
                    #if ( CG_SET[j]  ):
		    #CHARN[j] = n_groups
                    #	CG_SET[j] = 0

    # check unset heavy atom groups
    for pid_i, ptclObj_i  in strucC.ptclC:
        # If unset set chrge group
        #    if( CG_SET[i] ):
        if( ptclObj_i.tagsDict["qgroup"] == 0  ):
            NNAB = calc_nnab(pid_i,cov_nbindx)
            ELCNT = calc_elcnt(pid_i,strucC,cov_nblist,cov_nbindx)
            N_o = cov_nbindx[pid_i]
            N_f = cov_nbindx[pid_i+1] - 1
            if( ptclObj_i.tagsDict["number"] == 6 or ptclObj_i.tagsDict["number"] == 8 ):  #
		q_g = 0 
                for j_indx in range( N_o,N_f+1):
                    pid_j = cov_nblist[j_indx]
                    ptclObj_j = strucC.ptclC[pid_j]
                    if (ptclObj_j.tagsDict["qgroup"] != 0 ):
                        q_g = ptclObj_j.tagsDict["qgroup"]
                        break

		if( q_g  == 0 ):
			
		    n_groups = n_groups + 1
		    #CG_SET[i] = 0
		    ptclObj_i.tagsDict["qgroup"] = n_groups
		else:
		    #CG_SET[i] = 0
		    ptclObj_i.tagsDict["qgroup"] = q_g
		    
		#print " unset carbon set to ",ptclObj_i.tagsDict["qgroup"]
			

            if( ptclObj_i.tagsDict["fftype"] == 'CA' ):  # conjugate
                n_groups = n_groups + 1
                # CG_SET[i] = 0
                ptclObj_i.tagsDict["qgroup"] = n_groups
 
                    

    # set light atom groups 
    for pid_i, ptclObj_i  in strucC.ptclC:
        # If unset set chrge group
        #    if( CG_SET[i] ):
        if( ptclObj_i.tagsDict["qgroup"] == 0  ):
            NNAB = calc_nnab(pid_i,cov_nbindx)
            ELCNT = calc_elcnt(pid_i,strucC,cov_nblist,cov_nbindx)
            N_o = cov_nbindx[pid_i]
            N_f = cov_nbindx[pid_i+1] - 1
            
            if ( ptclObj_i.tagsDict["number"] == 1 or  ptclObj_i.tagsDict["number"] == 9 ):
                j_set = 0
                for j_indx in range( N_o,N_f+1):
                    pid_j = cov_nblist[j_indx]
                    ptclObj_j = strucC.ptclC[pid_j]
                    if (ptclObj_j.tagsDict["qgroup"] != 0 ):
                        n_groups = ptclObj_j.tagsDict["qgroup"]
                        j_set = 1
                if( j_set ):
                    ptclObj_i.tagsDict["qgroup"] = n_groups
                    #CG_SET[i] = 0
                else :
                    print ptclObj_i.tagsDict["symbol"] ,' connected to unset atom '
                    #print ptclObj_i.tagsDict["symbol"] , " - ",ptclObj_j.tagsDict["symbol"] 
                    #print ptclObj_i.tagsDict["fftype"] , " - ",ptclObj_j.tagsDict["fftype"] 
                    #for i in range(NA):
                    #    print i,ASYMB[i],ATYPE[i]
                    error_line = "!!!!!! Warning !!!!!!"
                    error_line += " unconnected atom "
                    print error_line
                    # sys.exit('charge group error')
                    
            if (ptclObj_i.tagsDict["fftype"] == 'LP' ):
                j_set = 0
                for j_indx in range( N_o,N_f+1):
                    pid_j = cov_nblist[j_indx]
                    ptclObj_j = strucC.ptclC[pid_j]
                    if (ptclObj_j.tagsDict["qgroup"] != 0 ):
                        n_groups = ptclObj_j.tagsDict["qgroup"]
                        j_set = 1
                if( j_set ):
                    ptclObj_i.tagsDict["qgroup"] = n_groups

                else :
                    print ' lone pair connected to unset atom '
                    print ptclObj_i.tagsDict["symbol"] , " - ",ptclObj_j.tagsDict["symbol"] 
                    print ptclObj_i.tagsDict["fftype"] , " - ",ptclObj_j.tagsDict["fftype"] 
                    sys.exit('charge group error')
                    

		
    debug = 0
    if( debug ):
        for pid_i, ptclObj_i  in strucC.ptclC:
            # If unset set chrge group
             print i,ASYMB[i],CG_SET[i],ptclObj_i.tagsDict["qgroup"],  RING_NUMB[i]
	     
	     
    # Check for unset atoms 
    for pid_i, ptclObj_i  in strucC.ptclC:
        if( ptclObj_i.tagsDict["qgroup"] == 0 ):
            NNAB = calc_nnab(pid_i,cov_nbindx)
            ELCNT = calc_elcnt(pid_i,strucC,cov_nblist,cov_nbindx)
            N_o = cov_nbindx[pid_i]
            N_f = cov_nbindx[pid_i+1] - 1
            
            print pid_i,ptclObj_i.tagsDict["symbol"],ptclObj_i.tagsDict["fftype"],NNAB,ELCNT[6],ELCNT[1],' not set '
            error_line = 'charge group error'
            print error_line
	    # sys.exit(error_line)
            #else :
            # print ASYMB[i],ATYPE[i],NNAB,ELCNT[6],ELCNT[1],GTYPE[i],' set '

    return 

def main():
    """
    Read in data file and replicate it 
    """

    #        
    # Read options 
    #
    options, args = get_options()
    #
    # Initialize mpi
    #
    p = mpiBase.getMPIObject()
    #
    # Read in fchk file into a simulation object 
    # 
    sim_o = SimulationGaussian(options.out_id)
    sim_o.readfchk(options.in_fchk)
    
    print sim_o
    write_xyz = True 
    if( write_xyz ):

        comment = " adf"
        append = False 
        sim_o.write_xmol("{}.xyz".format(options.out_id),comment,append)

    struc_o = sim_o.getstrucC()
    param_o = sim_o.getparamC()
    
    # initialize_fftags(struc_o)


    struc_o.build_bonded_nblist(max_nn=12.0,radii_buffer=1.25)
    struc_o.nblist_bonds()
    struc_o.nblist_angles()
    limdih = False 
    limitdih_n = 0
    struc_o.nblist_dih(limdih,limitdih_n)
    #
    # Use oplsaa types as fftype guess
    #   These should be checked and edited in the cply file
    #
    
    ring_nblist, ring_nbindex = find_rings(struc_o,struc_o.bonded_nblist, struc_o.bonded_nbindx)
    
    oplsaa_atomtypes(struc_o,struc_o.bonded_nblist, struc_o.bonded_nbindx)
    set_chargegroups(struc_o , ring_nblist, ring_nbindex,struc_o.bonded_nblist, struc_o.bonded_nbindx)

    print " Charge groups set "
    
    bb_o = Buildingblock()
    print " bb_o created  "
    bb_o.setStructureContainer(struc_o)
    print " bb_o strucC set to struc_o  "
    
    bb_o.set_cply_tags()
    # struc_o.zero_unitq()

    print " cply tags set    "
    
    bb_o.write_cply("{}.cply".format(options.out_id),write_ff=True,write_bonds=True)
    
    del sim_o


if __name__=="__main__":
    main()
