#! /usr/bin/env python
# Subroutines for reading and writing gaussian related files

# Dr. Travis Kemper
# NREL
# Initial Date 06/20/2013
# travis.kemper@nrel.gov

# Energy conversion
# http://physics.nist.gov/cgi-bin/cuu/Value?threv
HtoeV = 27.211385

def check_log(log_name):
    
    run_opt = 1
    try:
        with open(log_name) as f:
            print log_name, ' exist'
            F = open(log_name,'r')
            log_lines = F.readlines()
            F.close()
            # Check for a complete exicution
            if( len(log_lines) > 1):
                for line in log_lines:
                    last_line = line
                col = last_line.split()
                if ( len(col) > 1 ):
                    if ( col[0] == 'Normal' and col[1] == 'termination' ):
                        run_opt = 0
                        print '   and finished without error'
    except IOError:
       print ' gaussian output file ',log_name, ' does not exist needs to be generated '

    return run_opt

def read_optlog( log_name ):
    # Read log file of calculation using the opt keyword 

    debug = 0
    
    F = open(log_name,'r')
    Lines = F.readlines()
    F.close()
    #
    # Read in log file 
    #
    for line in Lines :
        col = line.split()
        if len(col) >= 7 :
            # Reset bonds for each print out 
            if ( col[1] == 'Name' and   col[2] == 'Definition' and   col[3] == 'Value' and   col[4] == 'Derivative' ):
                BONDS = []
                
                if(debug):
                    print " starting new bond array  "
            if col[2][:2] == 'R(' :
                b_id = col[2].split(',')
                b_id_i = b_id[0].split('(')
                b_i = int( b_id_i[1])
                b_id_j = b_id[1].split(')')
                b_j = int( b_id_j[0])
                BONDS.append([b_i-1,b_j-1])

                if(debug):
                    print b_i,b_j,b_i-1,b_j-1
    #
    # Read in log file 
    #
    for line in Lines :
        col = line.split()
        if len(col) >= 7 :
            # Reset bonds for each print out 
            if ( col[1] == 'Name' and   col[2] == 'Definition' and   col[3] == 'Value' and   col[4] == 'Derivative' ):
                ANGLES = []*3
            if col[2][:2] == 'A(' :
                id = col[2].split(',')
                id_i = id[0].split('(')
                i = int( id_i[1])
                j = int(id[1])
                id_k = id[2].split(')')
                k = int( id_k[0])                
                ANGLES.append([i-1,j-1,k-1])
    #
    # Read in log file 
    #
    for line in Lines :
        col = line.split()
        if len(col) >= 7 :
            # Reset bonds for each print out 
            if ( col[1] == 'Name' and   col[2] == 'Definition' and   col[3] == 'Value' and   col[4] == 'Derivative' ):
                DIH = []*4
            if col[2][:2] == 'D(' :
                id = col[2].split(',')
                id_i = id[0].split('(')
                i = int( id_i[1])
                j = int(id[1])
                k = int(id[2])
                id_l = id[3].split(')')
                l = int( id_l[0])
                DIH.append([i-1,j-1,k-1,l-1])
                
    
    return BONDS   , ANGLES, DIH 

def check_fchk( fchk_file ):
    
    bohr2angstrom = 0.5291772086

    run_calc = 1
    try:
        with open(fchk_file) as f:
            print fchk_file, ' exist'
            F = open(fchk_file,'r')
            Lines = F.readlines()
            F.close()
            # Check for a complete exicution
            if( len(Lines) > 1):
                for line in Lines:
				
		    col = line.split()
	    
		    if( len(col) > 2 ):
			if( col[0] == "Total" and col[1] == "Energy" ):
			    run_calc = 0
		    
		    
    except IOError:
       print ' gaussian output file ',fchk_file, ' does not exist needs to be generated '

	    
    return run_calc  


def parse_fchk( fchk_file ):
    import sys, numpy 
    
    bohr2angstrom = 0.5291772086

    R = []
    R_all = []
    ELN = []

    F = open(fchk_file,'r')
    Lines = F.readlines()
    F.close()

    read_r = 0
    read_eln = 0
    for line in Lines :
        col = line.split()

        if( read_r ):
            if (  col[0] == "Force" and col[1] == "Field" ):
                read_r = 0
                for atom_i in range(NA):
                    #print atom_i ,atom_i*3,atom_i*3+2,R_all[atom_i*3:atom_i*3+3]
		    vec_r_i =  numpy.array(  [R_all[atom_i*3:atom_i*3+3]] )
                    R.append(vec_r_i)
                
            else:
                for r_i in  map(float,col) :
                    R_all.append( r_i*bohr2angstrom )

        if( read_eln ):
            if ( len( ELN ) == NA ):
                read_eln = 0
            else:
                for eln_i in  map(int,col):
                    ELN.append( eln_i )

        
        if( len(col) > 2 ):
            if( col[0] == "Total" and col[1] == "Energy" ):
                TOTAL_ENERGY = float( col[3] )*HtoeV
            
        if( len(col) == 5 ):
            if( col[0] == "Number" and col[1] == "of"  and col[2] == "atoms" ):
                NA = int(col[4])
                
        if( len(col) == 6 ):
            if( col[0] == "Current" and col[1] == "cartesian"  and col[2] == "coordinates" ):
                read_r = 1
                
        if( len(col) > 2  ):
            if( col[0] == "Atomic" and col[1] == "numbers"   ):
                read_eln = 1
                
		
    return ( NA, ELN, R, TOTAL_ENERGY  )



def com_zmatrix(com_name):
    # Get zmatrix from gaussian input file

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

def  get_dih_id( zmatrix ):
    # find dihedral id's atoms and values from zmatrix 
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
                DIH_ID.append( colvar[0] )
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
		    DIH_ATOMS[dih_indx] =   [ zmat_l -1 , zmat_i-1 , zmat_j-1, zmat_k-1 ] 
                    
    return ( DIH_ID, DIH_VAL, DIH_ATOMS)
    
        
def print_com( id_name, ASYMB,R,ATYPE,CHARGES,ELECTRONS_i,options):
    from string import replace
    
    # Check multiplicity = 2 S + 1
    #   where S is the total spin of the system
    e_total = ELECTRONS_i - options.qm_charge 
    e_unpaired =  e_total % 2  
    if( e_unpaired == 0 ):
        multiplicity = 1 + options.qm_mult
    elif( e_unpaired == 1 ):
        multiplicity = 2 + options.qm_mult
        
        
    # set options 
    F_name = str(id_name) + ".com"
    chk_name = str(id_name) + ".chk"
    
    f = open(F_name,'w')
    f.write( "%s%s" % ("%chk=", chk_name))
    #f.write( "\n%snpros=%d " % ("%",options.qm_npros))
    f.write( "\n# P %s/%s  %s" % (options.qm_method,options.qm_basis,options.qm_kywd))
    f.write( "\n ")
    f.write( "\n %s " % (id_name))
    f.write( "\n ")
    f.write( "\n  %d %d " % (options.qm_charge,multiplicity))
    for i in range( len(ASYMB)):
        f.write("\n  %5s %16.8f %16.8f %16.8f " % (ASYMB[i],R[i][0],R[i][1],R[i][2]) )
    f.write( "\n ")
    f.write( "\n ")
    f.close()        
        


def com2zmat(com_id,comz_id,options):
    import sys, os
    from string import replace    
    #import shutil

    # Set .com name 
    com_name = str(com_id) + ".com"
    comz_name = str(comz_id) + ".com"
    
    # Copy current com to temp.com
    #shutil.copy2(comz_name, "temp.com")
    make_z = "newzmat -gencon  -redoz  "+ str( com_name)+ " "+str( comz_name)
    os.system(make_z)
    
    
    # modify com file to have correct chk file
    f = open(comz_name,'r')
    comz_lines = f.read()
    f.close()
    
    comz_lines = replace(comz_lines,com_id,comz_id)
    f = open(comz_name,'w')
    f.write(comz_lines)
    f.close
    
def run(options,calc_id):
    import sys, os
    from random import randint
    from string import replace 

    
    com_file = calc_id + ".com"
    log_file  = calc_id + ".log"
    random_numb = randint(2,9)


    # modify input file 
    f = open(com_file,'r')
    com_lines = f.read()
    f.close()

    if( options.cluster_host == "peregrine" ):
	load_gaussian = "module load gaussian/.g09_C.01 "
        user = 'tkemper'
        
            
        scratch_dir = "/gscratch1/"+ user+"/GAUSSIAN-"+calc_id+'-'+str(random_numb)
        mk_scratch = "mkdir " + scratch_dir
        rm_temp = "rm -rf /dev/shm/*"
        ex_scratch =  "export GAUSS_SCRDIR="+scratch_dir
        rm_scratch =  "rm -rf "+scratch_dir
    
        os.system(load_gaussian)
        os.system(mk_scratch)
        os.system(ex_scratch) 
        os.system(rm_temp)
        
        com_new = "%NoSave" 
        com_new = com_new + "\n" + "%RWF=/dev/shm/,1500MB,"+scratch_dir+"/,-1"
        com_new = com_new + "\n" + com_lines
        
        
    elif( options.cluster_host == "redmesa" ):
	load_gaussian = "module load gaussian/g09/C.01"
        user = 'twkempe'
        
        scratch_dir = "/gscratch1/"+ user+"/GAUSSIAN-"+calc_id+'-'+str(random_numb)
        mk_scratch = "mkdir " + scratch_dir
        rm_temp = "rm -rf /dev/shm/*"
        ex_scratch =  "export GAUSS_SCRDIR="+scratch_dir
        rm_scratch =  "rm -rf "+scratch_dir
    
        os.system(load_gaussian)
        os.system(mk_scratch)
        os.system(ex_scratch) 
        os.system(rm_temp)
            
        com_new = "%NoSave" 
        com_new = com_new + "\n" + "%RWF=/dev/shm/,1500MB,"+scratch_dir+"/,-1"
        com_new = com_new + "\n" + com_lines
    else:
        load_gaussian = ""
        com_new = com_lines

    mk_sm_subdir = 1
    work_dir = os.getcwd()                        
    if( mk_sm_subdir ):
        if (not os.path.isdir(calc_id)):
            os.mkdir(calc_id)
        os.chdir(calc_id)

        
    log_name = replace(com_file, "com", "log")

    # Print default route file
    f = file("Default.Route", "w")
    f.write("  -P-  %d " % (options.npros ))
    f.close()


    run_gaus=" g09 < "+str(com_file)+ " > "+ str(log_name)
    print run_gaus

    f = file(com_file, "w")
    f.write(com_new)
    f.close()

    os.system(run_gaus)
   
    debug = 0
    if( debug):
       os.system("pwd")
       sys.exit(' sp ')
       
    #Format chk
    if( check_log(log_name) == 0 ):
        chk2fchk( options , calc_id )

    os.chdir(work_dir)

    if( options.cluster_host == "peregrine" or options.cluster_host == "redmesa" ):
        os.system(rm_temp)
        os.system(rm_scratch)


        
def chk2fchk( options , calc_id ):
    import sys, os
    import file_io

    if( options.cluster_host == "peregrine" ):
	load_gaussian = "module load gaussian/.g09_C.01 "
        os.system(load_gaussian)
    elif( options.cluster_host == "redmesa" ):
	load_gaussian = "module load gaussian/g09/C.01"
        os.system(load_gaussian)
	
    chk_file = calc_id + ".chk"
    fchk_file = calc_id + ".fchk"

    chk_exists = file_io.file_exists(chk_file)
    if ( chk_exists ):
        make_fchk =    "formchk " + chk_file
        os.system(make_fchk)
        os.remove(chk_file) 
    else:
        fchk_exists = file_io.file_exists(fchk_file)
        if( fchk_exists ):
            print ' file ',fchk_file,' exists '
        else:
            print " no chk file ",chk_file
            sys.exit(' file not found ')
        