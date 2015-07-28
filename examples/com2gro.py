#! /usr/bin/env python

from structureContainer import StructureContainer
from particles     import Particle, ParticleContainer
import units

def get_options():
    import os, os.path
    from optparse import OptionParser
 
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    parser.add_option("--in_com", dest="in_com", type="string", default="", help="Input gaussain .com file ")
    parser.add_option("--out_gro", dest="out_gro", type="string", default="", help="Output .gro GROMACS file")
    #
    (options, args) = parser.parse_args()
        
    return options, args

def check_int(var):
    """
    check if variable var is an integer
    """
    
    try:
        x = int(var)
    except ValueError:
        return False
    
    return True
    
def read_com(strucC,data_file):
    """
    Read in structure information from gaussian com file

    Args:
        strucC (StructureContainer) 
        data_file (str) gaussian com file
    Return:
        strucC (StructureContainer) 
        
    """

    # Check to see if a previous read has occured
    pt_update = False
    #if( len(strucC.ptclC) > 0 ):
    #    pt_update = True

    F = open(data_file,'r')
    Lines = F.readlines()
    F.close()

    line_cnt = 0

    read_r = False
    
    for line in Lines :
        col = line.split()
        line_cnt += 1

        if( read_r ):
            if (  len(col) < 4  ):
                read_r = False
            else:
                p_i += 1
                symbol = str(col[0])
                r_i = [float(col[1]),float(col[2]),float(col[3]) ]
                
                if( pt_update ):
                    strucC.ptclC[p_i].position = r_i
                    strucC.ptclC[p_i].setTagsDict({"symbol":symbol})                    
                else:
                    pt_i = Particle( r_i ) 
                    pt_i.setTagsDict({"symbol":symbol})                    
                    strucC.ptclC.put(pt_i)

        if( len(col) == 2 ):
            if( check_int( col[0] ) and check_int( col[1] )  ):
                read_r = True
                p_i = 0
        
    return (strucC)

def write_gro(strucC,data_file):
    """
    Write gromacs structure file 
    """
    #
    latvec = strucC.getLatVec()

    gro_lines =  " com2gro \n"
    gro_lines += " %-2i \n"  %( int(len(strucC.ptclC)) )
    atom_indx = 0 
    for pid, pt_i  in strucC.ptclC:
        atom_indx += 1
        if( atom_indx > 10000): atom_indx = 1
        r_i = pt_i.position
        r_i_nm = [units.convert_angstroms_nm(r_i[0]) ,units.convert_angstroms_nm(r_i[1]) ,units.convert_angstroms_nm(r_i[2]) ]
        gro_lines += "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"  % ( pt_i.tagsDict["residue"],pt_i.tagsDict["resname"][:5],pt_i.tagsDict["gtype"][:5],atom_indx,r_i_nm[0],r_i_nm[1],r_i_nm[2] )
    gro_lines += " %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (units.convert_angstroms_nm(latvec[0][0]),units.convert_angstroms_nm(latvec[1][1]),units.convert_angstroms_nm(latvec[2][2]),units.convert_angstroms_nm(latvec[0][1]),units.convert_angstroms_nm(latvec[0][2]),units.convert_angstroms_nm(latvec[1][0]),units.convert_angstroms_nm(latvec[1][2]),units.convert_angstroms_nm(latvec[2][0]),units.convert_angstroms_nm(latvec[2][1])) 

    F = open( data_file, 'w' )
    F.write(gro_lines)
    F.close()

def add_gro_prop(strucC):
    """
     Add in particle properties need for .gro file
     """
    for pid, pt_i  in strucC.ptclC:
        add_dict = pt_i.tagsDict
        add_dict["residue"] = 1
        add_dict["resname"] = "MOLRES"
        add_dict["gtype"] = pt_i.tagsDict["symbol"]
        pt_i.setTagsDict(add_dict)                    
        
    return (strucC)

def main():
    """
    Read in gaussian fchk file and create an xyz file 
    """

    #        
    # Read options 
    #
    options, args = get_options()
    #
    #  Initialize blank system 
    # 
    struc_o = StructureContainer()
    #struc_o.unset_verbose()
    #param_o = ParameterContainer()

    struc_o = read_com(struc_o,options.in_com)
    struc_o = add_gro_prop(struc_o)
    write_gro(struc_o,options.out_gro)

    
if __name__=="__main__":
    main()

