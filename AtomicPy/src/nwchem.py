#! /usr/bin/env python
# IO for nwchem files 

# Dr. Travis Kemper
# NREL
# Initial Date 12/18/2013
# Email travis.kemper@nrel.gov
# Version 2.00

def write_et( ASYMB_i, R_i,Q_i, ASYMB_j, R_j,Q_j, out_et ):
    import sys 
    import elements , prop
    
    # ij n_i + c_j -> n_j + c_i
    
    ELN = []
    ELN.append( elements.asymb_eln(ASYMB_i) )
    ELN.append( elements.asymb_eln(ASYMB_j) )
    
    
    geom_ids = ["GEOMI","GEOMJ"]
    Q = [Q_i,Q_j]
    method = "UHF"
    
    f = open(out_et,"w")
    
    
    calc_id = "test"
    f.write( " start %s " % (calc_id))
    
    f.write( "\n geometry GEOMI units angstroms NOCENTER NOAUTOZ NOAUTOSYM" )
    for i in range( len(ASYMB_i)):
        f.write("\n  %5s %16.8f %16.8f %16.8f " % (ASYMB_i[i],R_i[i][0],R_i[i][1],R_i[i][2]) )
    f.write( "\n end" )
    
    f.write( "\n geometry GEOMJ units angstroms NOCENTER NOAUTOZ NOAUTOSYM" )
    for i in range( len(ASYMB_j)):
        f.write("\n  %5s %16.8f %16.8f %16.8f " % (ASYMB_j[i],R_j[i][0],R_j[i][1],R_j[i][2]) )
    f.write( "\n end" )

    # Print ij
    ASYMB_ij = ASYMB_i + ASYMB_j
    ELN_ij = ASYMB_i + ASYMB_j
    R_ij = R_i + R_j    
    f.write( "\n geometry GEOMIJ units angstroms NOCENTER NOAUTOZ NOAUTOSYM" )
    for i in range( len(ASYMB_ij)):
        f.write("\n  %5s %16.8f %16.8f %16.8f " % (ASYMB_ij[i],R_ij[i][0],R_ij[i][1],R_ij[i][2]) )
    f.write( "\n end" )

    f.write( "\n BASIS ")
    f.write( "\n * LIBRARY 6-31g ")
    f.write( "\n end ")

    for g_indx in range( len(geom_ids)):
        geom = geom_ids[g_indx]
        #print geom
        for q in Q:
            f.write( "\n SET geometry  %s " % (geom))
            f.write( "\n CHARGE  %d " % (q))
            f.write( "\n SCF ")
            f.write( "\n NOPEN 1 ")
            f.write( "\n UHF ")
            
            ELN_i = ELN[g_indx]
            M = prop.multiplicity(ELN_i,q)
            if( M == 1 ): MULTIPLICITY = "SINGLET"
            if( M == 2 ): MULTIPLICITY = "DOUBLET"
            if( M == 3 ): MULTIPLICITY = "TRIPLET"
            f.write( "\n %s " %  MULTIPLICITY )
            mofile = geom + "_"+str(q) + ".movecs"
            f.write( "\n VECTORS INPUT atom OUTPUT %s " %(mofile)) 
            
            f.write( "\n end ")
            f.write( "\n TASK SCF ")
    
    geom = "GEOMIJ"
    q = 1

    # Calculate i(n) + j(p)   
    q_i = 0
    molfile_i = geom_ids[0]  + "_"+str(q_i) + ".movecs"
    q_j = 1
    molfile_j = geom_ids[1]  + "_"+str(q_j) + ".movecs"
    molfile_ij = geom_ids[0]  + "_"+str(q_i) + geom_ids[1]  + "_"+str(q_j) + ".movecs"
    f.write( "\n SET geometry  %s " % (geom))
    f.write( "\n CHARGE  %d " % (q))
    f.write( "\n SCF ")
    f.write( "\n SYM OFF ")
    f.write( "\n ADAPT OFF ")
    f.write( "\n NOPEN 1 ")
    f.write( "\n %s " % (method))
    f.write( "\n VECTORS INPUT FRAGMENT %s %s OUTPUT %s " %(molfile_i,molfile_j,molfile_ij)) 
    f.write( "\n NOSCF ")    
    f.write( "\n end ")
    f.write( "\n TASK SCF ")
    

    # Calculate i(p) + j(n)   
    q_i = 1
    molfile_i = geom_ids[0]  + "_"+str(q_i) + ".movecs"
    q_j = 0
    molfile_j = geom_ids[1]  + "_"+str(q_j) + ".movecs"
    molfile_ji = geom_ids[0]  + "_"+str(q_i) + geom_ids[1]  + "_"+str(q_j) + ".movecs"
    f.write( "\n SET geometry  %s " % (geom))
    f.write( "\n CHARGE  %d " % (q))
    f.write( "\n SCF ")
    f.write( "\n SYM OFF ")
    f.write( "\n ADAPT OFF ")
    f.write( "\n NOPEN 1 ")
    f.write( "\n %s " % (method))
    f.write( "\n VECTORS INPUT FRAGMENT %s %s OUTPUT %s " %(molfile_i,molfile_j,molfile_ji)) 
    f.write( "\n NOSCF ")    
    f.write( "\n end ")
    f.write( "\n TASK SCF ")
    
    # ET calc
    f.write( "\n SET geometry  %s " % (geom))
    f.write( "\n CHARGE  %d " % (q))
    f.write( "\n ET ")
    f.write( "\n VECTORS REACTANTS %s " %(molfile_ij)) 
    f.write( "\n VECTORS PRODUCTS %s " %(molfile_ji)) 
    f.write( "\n end ")
    f.write( "\n TASK SCF ET")
    
    f.close()
    
    

def check_log(nw_log):
    import sys
    import file_io
    
    log_finished = 0 

    if ( file_io.file_exists(nw_log) ): 
        F = open(nw_log,'r')
        Lines = F.readlines()
        F.close()

        # Check to see if log file finished 
        for line in Lines :
            col = line.split()
            if( len(col) == 1 ):
                if( col[0].strip() == "CITATION" ):
                    log_finished = 1

    return      log_finished
    
def read_log(nw_log,geom_name):
    import sys ,string
    import file_io

    if ( file_io.file_exists(nw_log) ): 
        F = open(nw_log,'r')
        Lines = F.readlines()
        F.close()
    else:
        print nw_log
        sys.exit(' log file does not exist')

    # Check to see if log file finished 
    log_error = 1
    for line in Lines :
        col = line.split()
        if( len(col) == 1 ):
            if( col[0].strip() == "CITATION" ):
                log_error = 0

    if( log_error ):
        print nw_log
        sys.exit(' error in log file ')

        ASYMB = []
        R = []
    
    # Check to see if log file finished 
    log_error = 1
    line_cnt = 0
    read_geom = 0 
    for line in Lines :
        col = line.split()
        line_cnt += 1
        if( len(col) > 3 ):
            if(  col[0]  == "Geometry" and col[1].replace('"', '').strip() == geom_name and col[2] == "->"  ):
                read_geom = 1
                line_cnt = 0
            
                ASYMB = []
                R = []

        if( line_cnt >= 7 and read_geom and len(col) >= 5 ):
            ASYMB.append( col[1] )
            R.append( [float(col[3]),float(col[4]),float(col[5])] )

        if( read_geom and len(col) < 1 and line_cnt >= 7  ):
            read_geom = 0
            
        
    return ASYMB,R


def read_etlog(nw_log):
    import sys ,string
    import file_io

    if ( file_io.file_exists(nw_log) ): 
        F = open(nw_log,'r')
        Lines = F.readlines()
        F.close()
    else:
        print nw_log
        sys.exit(' log file does not exist')

    # Check to see if log file finished 
    log_error = 1
    for line in Lines :
        col = line.split()
        if( len(col) == 1 ):
            if( col[0].strip() == "CITATION" ):
                log_error = 0

    if( log_error ):
        print nw_log, " did not finish "
        sys.exit(' error in log file ')
    
    # Check to see if log file finished 
    log_error = 1
    read_et = 0
    vij = 0.0 
    for line in Lines :
        col = line.split()
        if( len(col) > 3 ):
            if(  col[0]  == "Electron" and  col[1]  == "Transfer"  and  col[2]  == "Coupling"   and  col[3]  == "Energy" ):
                read_et = 1

        if( read_et and len(col) > 0  and col[1] == "eV" ):
            vij = float(col[0])
            read_et = 0
            
    return vij
