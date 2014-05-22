#! /usr/bin/env python
"""
Scripts for hpc related opterations
"""

# Dr. Travis Kemper
# NREL
# Initial Date 10/24/2013
# travis.kemper@nrel.gov



def write_pbs(pbs_templ,calc_id,input_file,options):
    """
    Write pbs submission file 
    """
    from string import replace

    # Print pbs script
    pbs_id = calc_id+'.pbs'
    pbs_name = pbs_id
    pbs_dih = pbs_templ
    pbs_dih = replace(pbs_dih,"<calc_id>",calc_id)
    pbs_dih = replace(pbs_dih,"<input_file>",input_file)
    pbs_dih = replace(pbs_dih,"<pmem>",str(options.pmem) )
    pbs_dih = replace(pbs_dih,"<nnodes>",str(options.nnodes) )
    pbs_dih = replace(pbs_dih,"<npros>",str(options.npros) )
    f = file(pbs_name, "w")
    f.write(pbs_dih)
    f.close()
    
    return pbs_id




def submit_job( struct_dir, pbs_id ,options):
    """
    Submit pbs file to queue 
    """
    import sys, os

    print " submitting job " , pbs_id
    #submit = "sbatch " + pbs_id
    submit = options.submit_command +" "+ pbs_id

    #print " sumitting ",submit
    #sys.exit(' submit_job ')
    
    os.system(submit)
 
    return 
    #sys.exit('test job')
