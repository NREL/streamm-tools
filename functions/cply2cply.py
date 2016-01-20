import mpiBase
import sys 
from structureContainer import StructureContainer
from bonds         import Bond,     BondContainer
from buildingblocks import Buildingblock
from periodictable import periodictable

import numpy as np

def get_options():
    import os, os.path
    from optparse import OptionParser
 
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    parser.add_option("--in_cply", dest="in_cply", type="string", default="", help="Input .cply reference file")
    parser.add_option("-o","--out_id", dest="out_id", type="string", default="out", help="id")
    parser.add_option("-l","--setlabel", dest="setlabel", default=False,action="store_true", help="Re set particle labels ")
    parser.add_option("-m","--setmass", dest="setmass", default=False,action="store_true", help="Re set particle masses ")
    parser.add_option("-q","--setqgroup", dest="setqgroup", default=False,action="store_true", help="Find charge groups  ")
    parser.add_option("-r","--setrings", dest="setrings", default=False,action="store_true", help="Find rings ")
    parser.add_option("-s","--setresidues", dest="setresidues", default=False,action="store_true", help="Find residue numbers ")
    parser.add_option("-f","--setfftype", dest="setfftype", default=False,action="store_true", help="Find fftypes ")
    parser.add_option("-b","--rebond", dest="rebond", default=False,action="store_true", help="Rebond particles in cply file ")
    #
    # Filters
    #
    parser.add_option("--symbol", dest="symbol", type="string", default="", help=" select atoms of group by (atomic) symbol   ")    
    parser.add_option("--label", dest="label", type="string", default="", help="select atoms of group by label  ")    
    parser.add_option("--fftype", dest="fftype", type="string", default="", help="select atoms of group by force field type  ")    
    parser.add_option("--chain", dest="chain", type="string", default="", help="select atoms of group by chain/molecule number  ")    
    parser.add_option("--resname", dest="resname", type="string", default="", help="select atoms of group by residue name  ")    
    parser.add_option("--residue", dest="residue", type="string", default="", help="select atoms of group by resudue number  ")    
    parser.add_option("--ring", dest="ring", type="string", default="", help="select atoms of group by particlesn a ring   ")    
    
    (options, args) = parser.parse_args()
        
    return options, args


def addtagDic(dic_i,tag,tag_str,setint=False):
    """
    Take a string from input split it into values and add it to a dictionary list
    """
    if( len( tag_str ) ):
        dic_i[tag] = []
        for id_s in tag_str.split():
            if( setint ):
                dic_i[tag].append(int(id_s))
            else:
                dic_i[tag].append(id_s)
            
    return dic_i

def create_search(search_dic,f_symb,f_label,f_fftype,f_residue,f_resname,f_chain,f_ring):
    """
    Create a dictionary to pass to particle search
    """

    search_dic = addtagDic(search_dic,"symbol",f_symb)
    search_dic = addtagDic(search_dic,"label",f_label)
    search_dic = addtagDic(search_dic,"fftype",f_fftype)
    search_dic = addtagDic(search_dic,"residue",f_residue,setint=True)
    search_dic = addtagDic(search_dic,"resname",f_resname)
    search_dic = addtagDic(search_dic,"chain",f_chain,setint=True)
    search_dic = addtagDic(search_dic,"ring",f_ring,setint=True)
    
    return search_dic

def set_chargegroups(ff_charges,strucC , ring_nblist, ring_nbindex,cov_nblist, cov_nbindx ):
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
    pt = periodictable()

    #
    # Read in cply file into a structure container object
    # 
    bb_o = Buildingblock()
    bb_o.read_cply(options.in_cply)

    search_i = dict()
    search_i = create_search(search_i,options.symbol,options.label,options.fftype,options.residue,options.resname,options.chain,options.ring)
    if( len(search_i) > 0 ):
        if( options.verbose ):
            log_line = "\n Searching group i {} {} ".format(search_i,len(search_i))
            #log_out.write(log_line)
            if( options.verbose ): print log_line
        list_i = bb_o.ptclC.getParticlesWithTags(search_i)
        bb_i = bb_o.getSubStructure(list_i, particlesOnly=False )
    else:
        bb_i = bb_o
            
    if( len(bb_i.bondC) == 0 or options.rebond ):
        bb_i.bondC.clear()
        bb_i.ptclC.guess_radii()
        bb_i.build_bonded_nblist(max_nn=12.0,radii_buffer=1.25)
        bb_i.nblist_bonds()
    else:
        bb_i.bondC_nblist()

    if( options.setqgroup ): bb_i.set_chargegroups()
    bb_i.find_max_qgroup_id()
    if( options.setresidues ):  bb_i.set_residues()
    
    if( options.setlabel ): bb_i.set_label()
    if( options.setmass ): bb_i.set_mass()
    if( options.setrings ):  bb_i.find_rings()
    if( options.setfftype ):
        bb_i.oplsaa_atomtypes()
        set_biaryl = True 
        if(set_biaryl):
            bb_i.biaryl_types()
            bb_i.interring_types()

    bb_i.write_qgroup()
    
    # Set connecting termcap particles to lay along the x axis and the second connecter term to be at origin 
    bb_i.align_termcaps()
    print bb_i 
    comment = "Read in from {}  to output {} ".format(options.in_cply,options.out_id)
    append = False 
    bb_i.ptclC.write_xmol("{}.xyz".format(options.out_id),comment,append)
    bb_i.write_cply("{}.cply".format(options.out_id),write_ff=True,write_bonds=True)
    
    del bb_i


if __name__=="__main__":
    main()
