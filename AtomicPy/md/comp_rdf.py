# Compile multiple rdf files into a single file

def get_options():
    import os, os.path
    from optparse import OptionParser
 
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")

    parser.add_option("--rdf_out", dest="rdf_out", type="string", default="rdf.dat", help="Output rdf file ")

    (options, args) = parser.parse_args()
        
    return options, args
   
 
def main():
    #
    # Caclulate rdf
    #
    
    import os, sys, numpy , math 
    import datetime
    import  gromacs, elements, xmol, prop, file_io, groups  #, vectors 

    debug_read = 0  

    options, args = get_options()
    
    #
    # Loop over files and record inputs 
    #
    file_cnt = 0
    frame_cnt = 0
    box_den_j = []
    volume_i = [] #numpy.array()
    sum_i = 0
    sum_j = 0
    bin_size_flist = []
    r_cut_flist = []
    #
    for rdf_file in args:
	print "Getting inputs from ",rdf_file	
	if( file_io.file_exists(rdf_file)):
	    f_rdf = open(rdf_file,"r")
	    Lines = f_rdf.readlines()
	    f_rdf.close()
	    file_cnt += 1 
	    for line in Lines:
		col = line.split()
		if( len(col) >= 2 ):
		    if( col[0] == "#" and col[1] == "Frames" ):		    
			frame_cnt += int( col[2] )	
		    if( col[0] == "#" and col[1] == "N_i"  ):
			sum_i +=  int( col[2] ) 
		    if( col[0] == "#" and col[1] == "N_j"  ):
			sum_j +=  int( col[2] ) 
		    if( col[0] == "#" and col[1] == "Bin-size" ):		    
			bin_size_flist.append( float( col[2] ) )
		    if( col[0] == "#" and col[1] == "Cut-off" ):		    
			r_cut_flist.append( float( col[2] ) )
		if( len(col) >= 4 ):
		    if( col[0] == "#" and col[1] == "Box" and col[2] == "density"  and col[3] == "j"  ):
			box_den_j.append( float( col[4] ) )
		    if( col[0] == "#" and col[1] == "Average" and col[2] == "Box"  and col[3] == "Volume"  ):
			volume_i.append( float( col[4] ) )

    #
    # Allocate RDF array
    #
    bin_size = 0.0
    for sz_i in bin_size_flist :
	if ( sz_i > bin_size ): bin_size =  sz_i
    r_cut = 0.0
    for sz_i in r_cut_flist :
	if ( sz_i > r_cut ): r_cut =  sz_i
    #bin_size = numpy.maximum( bin_size_flist)      # Use largest 
    #r_cut = numpy.maximum( r_cut_flist)            # Use largest 
    n_bins = int(r_cut/bin_size)
    rdf_cnt = numpy.zeros(n_bins+1)
    #
    # Loop over frames and read in rdf counts 
    #
    for rdf_file in args:
	print " Reading counts from ",rdf_file	
	if( file_io.file_exists(rdf_file)):
	    f_rdf = open(rdf_file,"r")
	    Lines = f_rdf.readlines()
	    f_rdf.close()
	    for line in Lines:
		col = line.split()
		if( col[0] != "#" and  len(col) >= 3  ):
		    m_ij = float( col[1] )
		    bin_index = int( round( m_ij/bin_size) )
		    rdf_cnt[bin_index] += float( col[2] )
		    if( debug_read ):
			print m_ij,bin_index, float( col[2] )
		

    if( options.verbose ):
	print "      Finding averages "
    #
    # Find averages
    #
    box_vol_ave = numpy.average( volume_i )
    vol_cut = 4.0*math.pi/3.0*r_cut**3
    n_shperes = float(sum_i) #*float(file_cnt)
    total_cnts = numpy.sum( rdf_cnt)
    sphere_den_j = float(total_cnts)/vol_cut/n_shperes #/2.0  # N_B A^-3
    box_den_i = float(sum_i )/float(box_vol_ave)/float(file_cnt)
    box_den_j = float(sum_j )/float(box_vol_ave)/float(file_cnt)
    
    if( options.verbose ):
	print "   N_i ",sum_i
	print "   N_j ",sum_j
	print "   File_cnt ",file_cnt
	print "   Total counts ",total_cnts
	print "   Average box volume ",box_vol_ave
	print "   Volume of cut-off sphere ",vol_cut
	print "   Average box density i ",box_den_i," atoms/angstrom^3 "
	print "   Average box density j ",box_den_j," atoms/angstrom^3 "
	print "   Average cut-off sphere density ",sphere_den_j," atoms/angstrom^3 "
	
    # Write output 
    #
    rdf_file = options.rdf_out
    F_out = open(rdf_file,"w")
    F_out.write("# RDF Files  %s " %  (args))
    F_out.write("\n#    Bin-size %f  " % (bin_size))
    F_out.write("\n#    Cut-off %f  " % (r_cut))
    F_out.write("\n#    Files_cnt %d  " % (file_cnt))
    F_out.write("\n#    Total_cnts %d  " % (total_cnts))
    F_out.write("\n#    N_i %d " % (sum_i ))
    F_out.write("\n#    N_j %d " % (sum_j ))
    F_out.write("\n#    Average Box Volume %f " % ( box_vol_ave) )
    F_out.write("\n#    Box density i %f N A^-3 " % (box_den_i ))
    F_out.write("\n#    Box density j %f N A^-3 " % (box_den_j ))
    F_out.write("\n#    Sphere volume  %f A^3 " % (vol_cut ))
    F_out.write("\n#    Average Sphere density  %f N A^3 " % (sphere_den_j ))
    F_out.write("\n#    ")
    F_out.write("\n# bin index ; r     ; count_ave/frame ; dr vol ;  dr vol(aprox) ; g_sphere ; g_boxs  ")
    #                bin_index , r_val , dr_cnt_norm      , dr_vol,  dr_vol_apx,     sphere_g, box_g
    
    for bin_index in range( 1,n_bins):
	    
	r_val = bin_size*float(bin_index)
	r_in = r_val - bin_size*0.5
	r_out = r_val + bin_size*0.5

	dr_cnt = float( rdf_cnt[bin_index] )
	dr_cnt_norm =     dr_cnt    /float(file_cnt)
	dr_vol = 4.0*math.pi/3.0*( r_out**3 - r_in**3 )
	dr_vol_apx = 4.0*math.pi*(  r_val**2 )
	
	dr_rho = dr_cnt_norm/dr_vol
	sphere_g = dr_rho/sphere_den_j/float( sum_i )  *float(file_cnt)
	box_g = dr_rho/box_den_j/float( sum_i )  *float(file_cnt)
	
	F_out.write("\n  %d %f %f %f %f %f %f " % (bin_index,r_val,dr_cnt_norm,dr_vol,dr_vol_apx,sphere_g,box_g) )
	
    F_out.close()
	
	
if __name__=="__main__":
    main()
   