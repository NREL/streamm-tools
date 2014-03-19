#! /usr/bin/env python
# Read in geometry file and convert to topology

# Dr. Travis Kemper
# NREL
# Initial Date 05/07/2013
# travis.kemper@nrel.gov


#  Lists of atomic information for each atom in a structure  
#    ASYMB - list of atomic symbols
#    ELN   - list of atomic numbers
#    AMASS - list of atomic masses 


def get_options():
    import os, os.path
    from optparse import OptionParser
    usage = "usage: %prog <file to convert> [options] "
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    
    
    parser.add_option("--itp", dest="itp_file",  type="string", default="ff.itp",
                                    help="force field parameter file")
    parser.set_defaults(set_pmma=False)
    parser.add_option("--pmma", dest="set_pmma",action="store_true",
                                    help=" find pmma structures and set atom types accordingly ")
    parser.set_defaults(set_ptma=False)
    parser.add_option("--ptma", dest="set_ptma",action="store_true",
                                    help=" find ptma structures and set atom types accordingly ")
    parser.set_defaults(set_ions=False)
    parser.add_option("--ions", dest="set_ions",action="store_true",
                                    help=" find ions and set atom types accordingly ")
    
    parser.set_defaults(ff_charges=False)
    parser.add_option("--ff_charges", dest="ff_charges",action="store_true",
                                    help=" set charges to ff values described in atom_types.py ")

    parser.set_defaults(norm_dihparam=False)
    parser.add_option("--norm_dihparam", dest="norm_dihparam",action="store_true",help=" divide dihedral parameters in ff itp file by 4 ")
      
    parser.add_option("--exclusions", dest="exclusions",type="string", default="",
                                    help=" exclude some atom types, ids or symbols from read in ")
      
    parser.add_option("--sys_q", dest="sys_q",type=float, default=0.0,
                                    help=" System charge ")

    parser.add_option("--methyl_C", dest="methyl_C",type="string", default="C38",
                                    help="atom id of mythl carbon need for ptma ")
            
    parser.add_option("--out_gro", dest="out_gro", type="string", default="out.gro", help="gromacs output file ")
            
    (options, args) = parser.parse_args()

    return options, args


def main():
    #
    # Caclulate rdf
    #
    
    import os, sys, numpy , math 
    import datetime
    import  gromacs, top , pdb , elements, atom_types, lammps , xmol , gaussian
    
    debug_read = 0  

    options, args = get_options()

    ref_file = args[0]
    ff_file = options.itp_file

    LAT_CONST = numpy.zeros( (3,3) )
    
    LAT_CONST[0][0] = 50.0
    LAT_CONST[1][1] = 50.0
    LAT_CONST[2][2] = 50.0
    
    if( options.verbose ):
        print "Converting ", ref_file
        print "  using force field parameters found in ",ff_file

        
        # Print options
        if( options.set_ptma):
            print " ptma atypes will be located and set "
        else:
            print " no options have been specified "
            
    #
    # Get force field parameters
    #
    FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES = gromacs.read_itp(ff_file)
    

    const = []
    angle = []
    
    #
    # Define arrays
    #

    # Atomic properties 
    ELN = []
    CHARGES = []
    AMASS = []

    # structure  
    R = []*3
    BONDS = []*2
    ANGLES = []*3
    DIH = []*4
    IMPS = []

    #  FF data 
    ATYPE = []
    GTYPE = []
    RESID = []
    RESN = []
    CHARN = []

    # Ring info
    RINGLIST = []
    RINGINDEX  = []
    RING_NUMB = []
    
    #
    # Define reference file 
    ##
    #REF_FILE = log_file
    REF_SUF = ref_file[-4:]
    F = open(ref_file , 'r' )
    Lines = F.readlines()
    F.close()

    # New options that need to be passed 
    limdih =  0
    limitdih_n = 1
    
        
    if( REF_SUF == '.pdb'):        
        NA,ATYPE,GTYPE,RESID,RESN,CHARN,R,ASYMB,EXCLUSIONS,PDB_IND = pdb.atoms(Lines,options)

        ELN = elements.asymb_eln(ASYMB)
        AMASS = elements.eln_amass(ELN)
        RESID = top.initialize_resid( ELN )

        NBLIST, NBINDEX = pdb.nablist(Lines,EXCLUSIONS,PDB_IND,ATYPE)
        #NBLIST, NBINDEX = ff.build_nablist(ELN,R)

        debug = 0
        if( debug):
            for i in range( NA):                
                N_o = NBINDEX[ i  ]
                N_f = NBINDEX[  i+1  ] - 1
                print " atom ",i,i+1,ATYPE[i], N_f - N_o + 1
                for indx in range( N_o,N_f+1):
                    j = NBLIST[indx]
                    print "   connected to ",ATYPE[j],j
            sys.exit(" nab list debug")
        debug = 0
        
        BONDS = top.nblist_bonds(NA,NBLIST, NBINDEX)
        ANGLES = top.nblist_angles(NA,NBLIST, NBINDEX)
        DIH = top.nblist_dih(NA,NBLIST, NBINDEX,limdih,limitdih_n)
        IMPS = top.nblist_imp(NA,NBLIST, NBINDEX,ELN)
        
        CHARGES = top.initialize_charges( ELN )
        
    elif( REF_SUF == 'fchk'):

			    
	# Read in from fchk optimization 
	NA, ELN, R, TOTAL_ENERGY, CHARGES   = gaussian.parse_fchk( ref_file )
	ASYMB = elements.eln_asymb(ELN)

        AMASS = elements.eln_amass(ELN)
        RESID = top.initialize_resid( ELN )


        #   Build covalent nieghbor list for bonded information 
        NBLIST, NBINDEX = top.build_covnablist(ELN,R)
    
        BONDS = top.nblist_bonds(NA,NBLIST, NBINDEX)
        ANGLES = top.nblist_angles(NA,NBLIST, NBINDEX)
        DIH = top.nblist_dih(NA,NBLIST, NBINDEX,limdih,limitdih_n)
        IMPS = top.nblist_imp(NA,NBLIST, NBINDEX,ELN)
        
                
        #
        # Set GTYPES
        #
        ATYPE = []
        GTYPE = []
        RESN = []
        CHARN = []
        
        one = 1
        for i in range( len(ELN) ):
            CHARN.append(0)
            ATYPE.append(ASYMB[i])
            GTYPE.append(ASYMB[i])
            RESN.append(0)
    #
    # Set charge groups
    #
    CG_SET = []
    one = 1
    for i in range( len(ELN) ):
        CG_SET.append(one)
    #
    #
    d_mass = 0
    d_charge = 0
    RINGLIST, RINGINDEX , RING_NUMB = top.find_rings(ELN,NBLIST,NBINDEX)

    # Asign oplsaa atom types
    ATYPE, CHARGES = atom_types.oplsaa(  options.ff_charges,ELN,CHARGES,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB)
    
    #
    # Set charge groups
    #
    # CHARN = top.set_chargegroups(options,CG_SET,CHARN,ATYPE,ASYMB,ELN,R,NBLIST,NBINDEX, RING_NUMB,LAT_CONST)
    #clean_ids()
    
    # Optional atom type specificatons 
    if( options.set_ptma ):
        residue = 'PTMA'
        if( options.verbose ):
            print " setting PTMA "
        
        ATYPE,RESID,CHARGES,CG_SET,CHARN = atom_types.set_ptmatypes(  options ,ELN,ASYMB, ATYPE,GTYPE,RESID,CHARGES,AMASS,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB,CG_SET,CHARN )
        debug = 0
        if( debug):
            for i in range( len(IMPS) ):
                a_i = IMPS[i][0] + 1
                a_j = IMPS[i][1] + 1
                a_k = IMPS[i][2] + 1
                a_l = IMPS[i][3] + 1
                print ' dihedral index ',i, ATYPE[a_i-1],ATYPE[a_j-1],ATYPE[a_k-1],ATYPE[a_l-1],IMPTYPE_F[i]
            sys.exit('IMPTYPE_F')
            
            
                
  
    if( options.set_pmma ):
        residue = 'PMMA'
        ATYPE,RESID,CHARGES = atom_types.set_pmmatypes(options.sys_q, options.ff_charges, ELN, ATYPE,GTYPE,RESID,CHARGES,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB )
        
    if( options.set_ions ):
        residue = 'ION'
        print " setting ions "
        
        ATYPE,RESID,CHARGES = atom_types.set_iontypes(options.sys_q,  options.ff_charges,ELN, ATYPE,GTYPE,RESID,CHARGES,NBLIST,NBINDEX,RINGLIST, RINGINDEX , RING_NUMB )

    hack_rename = 1
    if( hack_rename ):

        if( options.verbose ):
            print "    Setting gromacs types to new index "
            
        h_cnt = 91
        lp_cnt = 0
        mol_atom = 0
        mol_cnt = 1
        atoms_mol = 500
        for i in range( len(ASYMB) ):
            a_symb = ASYMB[i].strip() 
            if ( GTYPE[i].strip() ==  a_symb):
                if( a_symb == 'H' ):
                    h_cnt = h_cnt + 1
                    GTYPE[i] = "H"+str(h_cnt)
                if( a_symb == 'LP' ):
                    lp_cnt = lp_cnt + 1
                    GTYPE[i] = "LP"+str(lp_cnt)
                print ' bad gtype ',[i]

    debug=0
    if( debug ):
        print NA
        for i in range(len(ELN)):            
            #print i+1,ASYMB[i] ,ELN[i],R[i][0],R[i][1],R[i][2],AMASS[i] ,ATYPE[i],CHARGES[i],GTYPE[i],RESID[i],RESN[i],CHARN[i]
            #print i+1 ,ELN[i],R[i][0],R[i][1],R[i][2],ATYPE[i],CHARGES[i],CHARGES[i],GTYPE[i],RESID[i],RESN[i],CHARN[i]
            print i+1 ,ELN[i],ASYMB[i],ATYPE[i],GTYPE[i],RESID[i],RESN[i],CHARN[i],CHARGES[i]
        sys.exit('print info')

        print ' bonds ', len(BONDS)
        for i in range( len(BONDS)):
            print i+1,BONDS[i][0]+1,BONDS[i][1]+1

        print ' angles ', len(ANGLES)
        for i in range( len(ANGLES)):            
            print ANGLES[i][0]+1,ANGLES[i][1]+1,ANGLES[i][2]+1
        print ' dihedrals ', len(DIH)                        
        for i in range( len(DIH)):
            print DIH[i][0]+1,DIH[i][1]+1,DIH[i][2]+1,DIH[i][3]+1
        
        sys.exit('print info')


    # Identify total number of atom types for lammps output 
    ATYPE_IND , ATYPE_REF,  ATYPE_MASS ,BTYPE_IND , BTYPE_REF, ANGTYPE_IND , ANGTYPE_REF, DTYPE_IND , DTYPE_REF = lammps.lmp_types(ELN,ATYPE,AMASS,BONDS,ANGLES,DIH)


    # Check atom types to be sure each atom of the same type has the same number of neighbors 
    #ATYPE_NNAB = top.check_types(ATYPE_IND , ATYPE_REF,GTYPE,NBLIST,NBINDEX)

    # Print output
    out_xyz = "out.xyz"
    xmol.write_xyz(ASYMB,R,out_xyz)
    xmol.print_cgxyz(CHARN,R)
    gromacs.print_gro(options.out_gro,GTYPE,RESID,RESN,R,LAT_CONST)
    
    
    ATYPE_NNAB = top.check_types(ATYPE_IND , ATYPE_REF,GTYPE,NBLIST,NBINDEX)    
    #ATYPE_EP, ATYPE_SIG = top.atom_parameters(ATYPE_IND , ATYPE_REF,  ATYPE_MASS,FF_ATOMTYPES)
    
    ATYPE_EP, ATYPE_SIG = top.atom_parameters(options.itp_file,ATYPE_IND , ATYPE_REF,  ATYPE_MASS,FF_ATOMTYPES)
    
    BONDTYPE_F , BONDTYPE_R0 ,BONDTYPE_K  = top.bond_parameters(options.itp_file,BTYPE_IND , BTYPE_REF,FF_BONDTYPES)
    ANGLETYPE_F , ANGLETYPE_R0 , ANGLETYPE_K = top.angle_parameters(options.itp_file,ANGTYPE_IND , ANGTYPE_REF,FF_ANGLETYPES)
    DIHTYPE_F ,DIHTYPE_PHASE ,DIHTYPE_K, DIHTYPE_PN,  DIHTYPE_C = top.dih_parameters(options.itp_file, options.norm_dihparam, DTYPE_IND , DTYPE_REF ,  FF_DIHTYPES,ATYPE_REF,ATYPE_NNAB  )
    
    IMPTYPE_F  = top.imp_parameters(options.itp_file)
    
    if( options.set_ptma ):        
        IMPS,IMPTYPE_F = atom_types.set_ptma_imps(NA,NBLIST, NBINDEX,ELN,ASYMB,IMPS,IMPTYPE_F)

    # Print new itp file with only used atom types and interactions
    new_itp = "ff-new.itp"
    
    norm_dihparam = 1
    # AT_LIST,NBD_LIST,ANG_LIST, DIH_LIST  = gromacs.print_itp(ASYMB,ATYPE,BONDS,ANGLES,DIH,IMPS,NBLIST,NBINDEX,FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES)
    AT_LIST,NBD_LIST,ANG_LIST, DIH_LIST  = gromacs.print_itp(new_itp,norm_dihparam,ASYMB,ATYPE,BONDS,ANGLES,DIH,IMPS,NBLIST,NBINDEX,FF_ATOMTYPES , FF_BONDTYPES , FF_ANGLETYPES ,  FF_DIHTYPES)
    #
    # set_chargeneutral()
    #
    #top_file = 'out_const.top'
    #gromacs.print_top( top_file, ASYMB , ELN,ATYPE, GTYPE, CHARN , CHARGES, AMASS,RESN, RESID ,BONDS , ANGLES , DIH , IMPS
    #                   ,DIH_CONST,DIH_CONST_ANGLE
    #                   ,BTYPE_IND, BONDTYPE_F, ANGTYPE_IND, ANGLETYPE_F
    #                   ,DTYPE_IND, DIHTYPE_F, IMPTYPE_F,LAT_CONST)    

    top_file = 'out.top'
    const = []
    angle = []
    gromacs.print_top( top_file, ASYMB , ELN,ATYPE, GTYPE, CHARN , CHARGES, AMASS,RESN, RESID ,BONDS , ANGLES , DIH , IMPS
                       ,const,angle
                       ,BTYPE_IND, BONDTYPE_F, ANGTYPE_IND, ANGLETYPE_F
                       ,DTYPE_IND, DIHTYPE_F, IMPTYPE_F,LAT_CONST)    

    #lammps.print_lmp(ATYPE_REF,ATYPE_MASS,ATYPE_EP,ATYPE_SIG,
    #          BTYPE_REF,BONDTYPE_R0,BONDTYPE_K,
    #          ANGTYPE_REF,ANGLETYPE_R0,ANGLETYPE_K,
    #          DIH,DTYPE_IND,DTYPE_REF,DIHTYPE_F,DIHTYPE_K,DIHTYPE_PN,DIHTYPE_PHASE,
    #          RESN,ATYPE_IND,CHARGES,R , ATYPE,
    #          BONDS ,BTYPE_IND, ANGLES ,ANGTYPE_IND, LAT_CONST)
    

    data_file = "out.data" 
    lammps.print_lmp(data_file,ATYPE_REF,ATYPE_MASS,ATYPE_EP,ATYPE_SIG,
          BTYPE_REF,BONDTYPE_R0,BONDTYPE_K,
          ANGTYPE_REF,ANGLETYPE_R0,ANGLETYPE_K,
          DIH,DTYPE_IND,DTYPE_REF,DIHTYPE_F,DIHTYPE_K,DIHTYPE_PN,DIHTYPE_PHASE,DIHTYPE_C,
          RESN,ATYPE_IND,CHARGES,R , ATYPE,
          BONDS ,BTYPE_IND, ANGLES ,ANGTYPE_IND, LAT_CONST)

    
if __name__=="__main__":
    main()
   
       