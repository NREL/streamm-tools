# Get data from remote cluster 


def get_options():
    import os, os.path
    from optparse import OptionParser
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)

    parser.add_option("-v","--verbose", dest="verbose", default=True, help="Verbose output ")
    
    # Cluster options
    parser.add_option("--cluster_host", dest="cluster_host",type="string",default="peregrine",help=" name of cluster ")

    # Getting data from cluster options
    parser.add_option("--info_dir", dest="info_dir",type="string",default="/hom/",help=" Directory on cluster where reference files are located ")
    parser.add_option("--info_file", dest="info_file",type="string",default="info_files.dat",help=" reference file name ")
    parser.add_option("--update_info", dest="update_info",type="string",default="True",help=" download reference file from cluster ")

    # should be reference file 
    parser.add_option("--qm_sufix", dest="qm_sufix",type="string",default="_qm2",help=" sufix of qm data file  ")
    
    #parser.add_option("--update_data", dest="update_data",type="int",default=1,help=" download data files specified in ")
    #parser.add_option("--update_plots", dest="update_plots",type="int",default=1,help=" update_plots ")
    #parser.add_option("--draw_str", dest="draw_str",type="int",default=1,help=" draw_str ")
    #parser.add_option("--update_vmd", dest="update_vmd",type="int",default=1,help=" update_vmd ")
    #parser.add_option("--update_gnuplot", dest="update_gnuplot",type="int",default=1,help=" update_gnuplot ")
    #parser.add_option("--open_files", dest="open_files",type="int",default=0,help=" open files")

    (options, args) = parser.parse_args()
    
    return options, args
    
def get_info(cluster_id,user_id,options):
    import sys, os
    
    get_dat = "scp "+ user_id+"@"+cluster_id + options.info_dir +"/"+ options.info_file +' ./ \n'
    print get_dat
    os.system(get_dat)


def get_data(cluster_id,user_id,struct_dir,job_name,dih_id,options):
    import sys, os
    
    debug = 0
    
    mk_dir = 'mkdir -p '+struct_dir+'/'
    os.system(mk_dir)
    
    dih_qm = struct_dir +'/' +job_name + '-' + dih_id + options.qm_sufix + '.dat'
        
    cluster_file = "scp "+ user_id+"@"+cluster_id + options.info_dir+"/" + dih_qm +' '+struct_dir+' \n' 
    get_dat = cluster_file
    
    if(debug):
        print get_dat        
        #sys.exit(" get_data ")
    else:
        os.system(get_dat)
    
def main():
    import sys, os 
    import string
    import file_io
    
    options, args = get_options()

    if( options.cluster_host == "peregrine" ):
        cluster_id = 'peregrine-login1.nrel.gov:'
        user_id = 'tkemper'
        
    elif( options.cluster_host == "redmesa" ):
        cluster_id = 'redmesa-login1.sandia.gov:'
        user_id = 'twkempe'
    else:
        print " unknow host id "
        sys.exit(" option error ")
    
        
    print " Downloading data files for ",user_id,"@",cluster_id


    if( options.update_info == "True" ):
        get_info(cluster_id,user_id,options)

    f = open(options.info_file,'r')
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
                
            print "   Getting ",struct_dir,job_name
            

            if( options.verbose ):
                print " Getting data from cluster " ,user_id,user_id,struct_dir,job_name,dih_id
            get_data(cluster_id,user_id,struct_dir,job_name,dih_id,options)

        
        if( len(col) >= 4 and col[0] == "ff_dih" ):
            
            struct_dir = col[1].strip()
            job_name = col[2].strip()
            ff_software = col[3].strip()
            ff_type_id = col[4].strip()

            dih_id = col[5].strip() 
            a_k  = int(  col[6].strip() )
            a_i  = int(  col[7].strip() )
            a_j  = int(  col[8].strip() )
            a_l  = int(  col[9].strip() )
            cent_min  = int(  col[10].strip() )
            cent_max = int(  col[11].strip() )
            cent_step = int(  col[12].strip() )
            
            print "   Getting ff data ",struct_dir,job_name
            

            if( options.verbose ):
                print " Getting data from cluster " ,user_id,user_id,struct_dir,job_name,dih_id
            get_data(cluster_id,user_id,struct_dir,job_name,dih_id,options)
        
            debug = 0
            
            mk_dir = 'mkdir -p '+struct_dir+'/'
            os.system(mk_dir)
            
            dih_ff = struct_dir +'/' +job_name + '-' + dih_id + "_ff" + ff_software + ff_type_id +".dat"

            cluster_file = "scp "+ user_id+"@"+cluster_id + options.info_dir +"/"+ dih_ff +' '+struct_dir+' \n' 
            get_dat = cluster_file
            
            if(debug):
                print get_dat        
                #sys.exit(" get_data ")
            else:
                os.system(get_dat)

            xmol_ff = struct_dir +'/' +job_name + '-' + dih_id + "_ff" + ff_software + ff_type_id +".xmol"

            cluster_file = "scp "+ user_id+"@"+cluster_id + options.info_dir +"/"+ xmol_ff +' '+struct_dir+' \n' 
            get_dat = cluster_file
            
            if(debug):
                print get_dat        
                #sys.exit(" get_data ")
            else:
                os.system(get_dat)



if __name__=="__main__":
    main()

