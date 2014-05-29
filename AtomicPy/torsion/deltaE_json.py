# Plot torsional data based on reference file 


def get_options():
    import os, os.path
    from optparse import OptionParser
    usage = "usage: %prog [options] [input_files] \n"
    parser = OptionParser(usage=usage)

    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")

    (options, args) = parser.parse_args()
    
    return options, args

def main():
    import sys, os , numpy 
    import string  
    import file_io, jsonapy
    from string import replace
    
    options, args = get_options()
    
    

    # Get lines of reference data file
    ref_data_file = args[0]
    if( file_io.file_exists(ref_data_file) ):
	    
	if( options.verbose ):
	    print " Reading in file ",ref_data_file
	f = open(ref_data_file,'r')
	ref_lines = f.readlines()
	f.close()
	
	ref_angle = []
	ref_energy = []
	
	for line in ref_lines:
	    col = line.split()
	    if( len(col) >= 4 and col[0] != "#" ):
		ref_angle.append( float( col[1] ) )
		ref_energy.append( float( col[2] ) )
    
	ref_min_energy = min(ref_energy)
	print " min ref energy ",ref_min_energy
	
	# Get lines of  data file
	data_file = args[1]
	if( file_io.file_exists(data_file) ):
		
	    if( options.verbose ):
		print " Reading in file ",data_file
	    f = open(data_file,'r')
	    lines = f.readlines()
	    f.close()
	    
	    angle = []
	    energy = []
	    
	    for line in lines:
		col = line.split()
		if( len(col) >= 4 and col[0] != "#" ):
		    angle.append( float( col[1] ) )
		    energy.append( float( col[2] ) )
	
	    min_energy = min(energy)
	    print " min  energy ",min_energy
	    
	    dEdata_file = args[2]
	    
	    dE_file = open( dEdata_file, 'w')
	    dE_file.write("# dE = %s - %s " % (data_file,ref_data_file))
	    dE_file.write("\n# file 2  %s " % (data_file))
	    dE_file.write("\n# file 1 %s " % (ref_data_file))
	    dE_file.write("\n# angle (deg), energy file 2 (eV), energy file 1 (eV), delta E (eV) " )
	    # loop over data and compare energy to reference 
	    for indx in range( len(angle) ):
		angle_i = angle[indx]
		    
		# loop over reference to find same angle
		for ref_indx in range( len(ref_angle) ):
		    ref_angle_i = ref_angle[ref_indx]
		    if( ref_angle_i == angle_i ):
			delta_energy_i = energy[indx] - min_energy - ref_energy[ref_indx]  + ref_min_energy
			#print  energy[indx] ,   ref_energy[ref_indx] ,delta_energy_i
			dE_file.write("\n %f8 %f16 %f16 %f16" % (angle_i,energy[indx] - min_energy , ref_energy[ref_indx]  - ref_min_energy,delta_energy_i))
			
	    dE_file.close()
	    
	    # Append gnuplot file
	    calc_i = 1
	    plot_l =  " \'"+dEdata_file+"\' " + ' us 1:($4)*EVKC  axis x1y1 w l ls '+str(calc_i+1)+" title \'n=\'  smooth unique  , \\" + "\n"
	    gnuplot_file = 'dE_D100.plot'
	    gnuplot_f = open(gnuplot_file,'a')
	    gnuplot_f.write(plot_l)
	    gnuplot_f.close()
	    
if __name__=="__main__":
    main()

