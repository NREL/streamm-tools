"""
Class of elements with associated atomic properties 
"""

# Dr. Travis Kemper
# NREL
# Initial Date 7/01/2014
# travis.kemper@nrel.gov



class element():
    """
    Gives the basic properties of the elements in the periodic table
    """


    def __init__(self,symbol,number,mass,ionic_radii):
        """
        Constructor for a element
        
        Args:
        symbol,number,ionic_radii
        """

        if isinstance(symbol, str):
            self.symbol = symbol
        else:
            print "1st arg should be str"
            raise TypeError


        if isinstance(number, int):
            self.number = number
        else:
            print "1st arg should be int"
            raise TypeError

        if isinstance(mass, float):
            self.mass = mass
        else:
            print "1st arg should be float"
            raise TypeError


        if isinstance(ionic_radii, float):
            self.ionic_radii = ionic_radii
        else:
            print "1st arg should be float"
            raise TypeError


    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.symbol 
        del self.number 
        del self.mass 
        del self.ionic_radii

    def symbol(self):
        return self.symbol


    def mass(self):
        return self.mass

class periodictable:
    """
    Elements of the periodic table 
    """

    def __init__(self):
        """
        Set elements of the periodic table 
        """
        
        self.elements=dict()
        self.maxgid = 110

        symbol_list = "H"
        number_list = 1
        ionic_radii = 1.01

        H = element("H",1,1.000,1.01)

        self.elements["H"] = H

        C = element("C",6,12.00,2.01)

        self.elements["C"] = C


    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.elements
        del self.maxgid

    def getelementWithMass(self,mass_i):
        """
        Find element based on mass
        """


        mass_i_int = int(mass_i)
        for el_symb in self.elements:
            el_mass_int = int( self.elements[el_symb].mass )
            if( mass_i_int == el_mass_int ):
                el =  self.elements[el_symb]

        return el

    
