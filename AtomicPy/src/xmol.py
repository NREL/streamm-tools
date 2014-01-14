#! /usr/bin/env python
# subroutines for reading and writing xmol style files

# Dr. Travis Kemper
# NREL
# Initial Date 10/09/2013
# travis.kemper@nrel.gov

def write_xyz(ASYMB,R,out_xyz):
    # Write file in xyz or single frame xmol file format 
    
    import sys 
    #
    # printing xmol 
    #
    F = open( out_xyz, 'w' )
    F.write( "%d " % ( len(ASYMB) )  )
    F.write( '\n' )
    for i in range( len(ASYMB) ):
        F.write("\n  %5s %16.8f %16.8f %16.8f " % (ASYMB[i],R[i][0],R[i][1],R[i][2]) )
    F.close()

    
def print_cgxyz(CHARN,R):
    import sys 
    print
    print ' printing charge group xyz '
    print
    spc = str(' ')
    xyz_form = "%5s %16.8f %16.8f %16.8 "
    F = open( 'chargegroup.xyz' , 'w' )
    F.write( "%d \n" % ( len(CHARN) )  )
    F.write( '\n' )
    for i in range( len(CHARN) ):
        i = int(i)
        # print i+1,ASYMB[i],R[i][0],R[i][1],R[i][2]
        # F.write( xyz_form % (ASYMB[i],R[i][0],R[i][1],R[i][2] ) )
        chr_n = CHARN[i]
        if( chr_n > 100 ): chr_n = chr_n - int( chr_n/100 ) * 100
        line = str(chr_n ) +spc+ str(R[i][0])+spc+str(R[i][1])+spc+str(R[i][2])
        F.write( line + '\n' )
    F.close()


def print_xmol(ASYMB,R,file_xmol):
    import sys , os
    import file_io
    #
    debug = 0
    #
    # appending xmol
    #
    current_frame =  str( "%d \n" % ( len(ASYMB) )  )
    current_frame = current_frame + str( '\n' )
    for i in range( len(ASYMB) ):
        current_frame = current_frame + str( "%5s %16.8f %16.8f %16.8f \n"  % (ASYMB[i] ,R[i][0], R[i][1],R[i][2]) )

    if( file_io.file_exists( file_xmol ) ):

        if( debug): print " File exists and will be appended "
        str_file = open( file_xmol, 'r' )
        prev_frames = str_file.read()
        str_file.close()
        all_frames = prev_frames + current_frame
    else:
        if( debug): print " File does not exists "
        all_frames = current_frame

    F = open( file_xmol, 'w' )
    F.write( all_frames )
    F.close()

