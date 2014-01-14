#! /usr/bin/env python
# run ab initio torsional potential 

# Dr. Travis Kemper
# NREL
# 12/09/2013
# travis.kemper@nrel.gov



def main():
    import sys, os , string 
    import file_io, gaussian
    
    
    
    # Store working dir  
    work_dir = os.getcwd()
    
    # Read index files from args
    for indx_file in args:
        # Get lines of index file   
        f = open(indx_file,'r')
        Lines = f.readlines()
        f.close()
        for line in Lines:
            col = line.split()
            if( len(col) >= 4 and col[0] == "qm_dih" ):
                
                mol_dir = col[1].strip()
                mol_id = col[2].strip()
                mol_repeat = int(col[3])
                mol_acc = col[4].strip()
                
                # File info
                struct_dir = mol_dir + "/" + mol_id + "/"
                job_name = mol_acc + "_" + mol_id + "_n" + str(mol_repeat)
                
                
                
		dih_id = col[5].strip()
		#a_k,a_i, a_j ,a_l,
		cent_min  = int(  col[10].strip() )
		cent_max = int(  col[11].strip() )
		cent_step = int(  col[12].strip() )
                
                
                print struct_dir , job_name , dih_id , cent_min , cent_max , cent_step