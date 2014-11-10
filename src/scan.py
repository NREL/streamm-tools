"""
Potential energy scan object
"""

import numpy as np

class PES():
    """
    Potential energy surface
    """

    def __init__(self,method,coord_list,energy_list,method_list,basis_list,tag_list, verbose=False):
        """
        Constructor: for pes class

        Args:
            method (str) method used to calculate pes 
            verbose (bool): flag for printing status/debug info
        """
        self.verbose=verbose        
        self.method = method
        self.coord_list = np.array(coord_list)        # Coordinate point of energy calculation 
        self.energy_list = np.array(energy_list)    # Energy at coordinate 
        self.method_list = method_list
        self.basis_list = basis_list
        self.tag_list = tag_list

        self.min_energy_list = self.calc_min_energy_list()     # shifted to min
        self.first_energy_list = self.calc_first_energy_list()      # shifted to zero coord

    def calc_min_energy_list(self):
        """
        Calculate the pes shifted above zero
        """
        return self.energy_list - self.min_energy()

    def calc_first_energy_list(self):
        """
        Calculate the pes shifted to first point 
        """
        return self.energy_list - self.energy_list[0]
        
    def __del__(self):
        """
        Destructor, clears dictionary memory
        """
        del self.method
        del self.energy_list
        del self.coord_list
        del self.min_energy_list
        del self.first_energy_list 
        del self.method_list 
        del self.basis_list 
        del self.tag_list
        
    def __len__(self):
        """
        'Magic' method for returning size of container
        """
        return len(self.energy_list)

    def min_energy(self):
        """
        Find minimum energy
        """
        return min(self.energy_list)
    

    def __str__(self):
        """
        'Magic' method for printng contents of container
        """
        str_line = " # coordinate ; energy \n"
        for indx in range(len(self.energy_list)):
            str_line += " %f %f  %f %f \n"%(self.coord_list[indx],self.energy_list[indx],self.min_energy_list[indx],self.first_energy_list[indx] )
        return str_line

    def get_coord(self):
        """
        Return coord list
        """
        return self.coord_list

    def get_energy(self):
        """
        Return energy list
        """
        return self.energy_list


    def get_min_energy(self):
        """
        Return min energy list
        """
        return self.min_energy_list

    def write(self,out_pes):
        """
        Write file of pes
        """
        
        str_line = " # coordinate ; energy ; method ; basis ; tag  \n"
        for indx in range(len(self.energy_list)):
            str_line += "%d %f %f  %s %s %s \n"%(indx,self.coord_list[indx],self.energy_list[indx],self.method_list[indx],self.basis_list[indx],self.tag_list[indx] )
        F = open(out_pes,"w")
        F.write(str_line)
        F.close()


def delta(pes_1,pes_2):
    """
    Take the difference in energy between two pes
    """
    verbose = True

    # Initialize PES lists 
    coord_list = []
    energy_list = []
    method_list = []
    basis_list = []
    tag_list = []

    for indx in range(len(pes_1.energy_list)):
        if( pes_1.coord_list[indx] != pes_2.coord_list[indx] ):
            error_line = " Coordinate %d of PES's %f and %f don't match "%( indx,pes_1.coord_list[indx],pes_2.coord_list[indx] )
            sys.exit(error_line)

        if( verbose ):
            print " %f %f  %f %f \n"%(pes_1.coord_list[indx],pes_1.min_energy_list[indx],pes_2.coord_list[indx],pes_2.min_energy_list[indx] )

        coord_list.append(pes_1.coord_list[indx])
        energy_list.append(pes_1.min_energy_list[indx] - pes_2.min_energy_list[indx])
        method_list.append( "%s-%s"%(pes_1.method_list[indx] ,pes_2.method_list[indx]) )
        basis_list.append("%s-%s"%(pes_1.basis_list[indx] ,pes_2.basis_list[indx]) )
        tag_list.append("%s-%s"%(pes_1.tag_list[indx] ,pes_2.tag_list[indx]) )

    return PES("delta",coord_list,energy_list,method_list,basis_list,tag_list )



def add(pes_1,pes_2):
    """
    Add two PES 
    """
    verbose = True

    # Initialize PES lists 
    coord_list = []
    energy_list = []
    method_list = []
    basis_list = []
    tag_list = []

    for indx in range(len(pes_1.energy_list)):
        if( pes_1.coord_list[indx] != pes_2.coord_list[indx] ):
            error_line = " Coordinate %d of PES's %f and %f don't match "%( indx,pes_1.coord_list[indx],pes_2.coord_list[indx] )
            sys.exit(error_line)

        if( verbose ):
            print " %f %f  %f %f \n"%(pes_1.coord_list[indx],pes_1.min_energy_list[indx],pes_2.coord_list[indx],pes_2.min_energy_list[indx] )

        coord_list.append(pes_1.coord_list[indx])
        energy_list.append(pes_1.min_energy_list[indx] + pes_2.min_energy_list[indx])
        method_list.append( "%s+%s"%(pes_1.method_list[indx] ,pes_2.method_list[indx]) )
        basis_list.append("%s+%s"%(pes_1.basis_list[indx] ,pes_2.basis_list[indx]) )
        tag_list.append("%s+%s"%(pes_1.tag_list[indx] ,pes_2.tag_list[indx]) )

    return PES("delta",coord_list,energy_list,method_list,basis_list,tag_list )


def read(in_pes):

    """
    Read file of pes
    """
    verbose = True 

    coord_list = []
    energy_list = []
    method_list = []
    basis_list = []
    tag_list = []

    F = open(in_pes , 'r' )
    en_lines = F.readlines()
    F.close()
    if( verbose ): print " %s energy file found "%(in_pes)
    for en_line in en_lines:
        col = en_line.split()
        if( col[0] != "#" and len(col) >= 4 ):
            coord_list.append(float(col[1]))
            energy_list.append( float(col[2]))
            method_list.append(col[3])
            basis_list.append(col[4])

            if( len(col) >= 6 ):
                tag = str(col[5:])
            else:
                tag = "NA"

            tag_list.append( tag)

    method = method_list[0]
    
    return PES(method,coord_list,energy_list,method_list,basis_list,tag_list )
