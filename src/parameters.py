"""
Class data structures force-field parameters 
"""
import copy, sys

class ljtype:
    """
    Set of LJ parameters 
    """

    def __init__(self, ptype1="blank"):
        """
        Constructor for a Lennard-Jones parameter.
        
        Args:
             ptype1  (str)   Atom type 
        """

        if isinstance(ptype1, str):
            self.ptype1 = ptype1
        else:
            print "1st arg should be str"
            raise TypeError


        # Set default values for parameters
        self.epsilon = 0.0
        self.sigma = 0.0 

    def __del__(self):
        """
        Destructor, clears object memory
        """

        
    def __str__(self):
        """
        'Magic' method for printng contents
        Delegates to __str__ method for contained objects
        """
        strucStr =  " LJ  %s  epsilon %f sigma %f "%(self.ptype1, self.epsilon, self.sigma )
                
        return strucStr

    def setparam(epsilon,sigma):
        """
        set Harmonic parameters

        E = 4 epsilon [ (sigma/r_ij)^12 - (sigma/r_ij)^6]

        Args:
            epsilon (float)  energy  kcal/mol  
            sigma    (float) distance
        """
        

        if isinstance(epsilon, float):
            self.epsilon = epsilon
        else:
            print "1st arg should be float"
            raise TypeError

        if isinstance(sigma, float):
            self.sigma = sigma
        else:
            print "2nd arg should be float"
            raise TypeError

        

class bondtype:
    """
    Set of bond parameters 
    """

    def __init__(self, ptype1="blank", ptype2="blank" , type="harmonic" ):
        """
        Constructor for a bond parameter.
        
        Args:
             ptype1  (str)   Atom type 
             ptype2  (str)   Atom type 
             type    (str)   Bond type 
        """

        if isinstance(ptype1, str):
            self.ptype1 = ptype1
        else:
            print "1st arg should be str"
            raise TypeError

        if isinstance(ptype2, str):
            self.ptype2 = ptype2
        else:
            print "2nd arg should be str"
            raise TypeError

        if isinstance(type, str):
            self.type = type
        else:
            print "3rd arg should be str"
            raise TypeError

        # Set default values for parameters
        self.r0 = 0.0
        self.kb = 0.0 

    def __del__(self):
        """
        Destructor, clears object memory
        """

        
    def __str__(self):
        """
        'Magic' method for printng contents
        Delegates to __str__ method for contained objects
        """
        strucStr =  " bond  %s - %s type %s "%(self.ptype1,self.ptype2,self.type)
        
        if( self.type ==  "harmonic" ):
            strucStr += "\n  harmonic r_0 = %f K = %f " %(self.r0 ,self.kb )
                
        return strucStr

    def setharmonic(self, r0, kb):
        """
        set Harmonic parameters

        E = kb( r - r0 )^2 

        Args:
            r0 (float) distance 
            kb (float) force constant 
        """

        
        if isinstance(r0, float):
            self.r0 = r0
        else:
            print "1st arg should be float"
            raise TypeError

        if isinstance(kb, float):
            self.kb = kb
        else:
            print "2nd arg should be float"
            raise TypeError
    
        
class angletype:
    """
    Set of Angle parameters 
    """

    def __init__(self, ptype1="blank", ptype2="blank", ptype3="blank" , type="harmonic" ):
        """
        Constructor for a angle parameter.
        
        Args:
             ptype1  (str)   Atom type 
             ptype2  (str)   Atom type 
             ptype3  (str)   Atom type 
             type    (str)   Bond type 
        """

        if isinstance(ptype1, str):
            self.ptype1 = ptype1
        else:
            print "1st arg should be str"
            raise TypeError

        if isinstance(ptype2, str):
            self.ptype2 = ptype2
        else:
            print "2nd arg should be str"
            raise TypeError

        if isinstance(ptype3, str):
            self.ptype3 = ptype3
        else:
            print "3rd arg should be str"
            raise TypeError

        if isinstance(type, str):
            self.type = type
        else:
            print "4th arg should be str"
            raise TypeError

        # Set default values for parameters
        self.theta0 = 0.0
        self.kb = 0.0 

    def __del__(self):
        """
        Destructor, clears object memory
        """

        
    def __str__(self):
        """
        'Magic' method for printng contents
        Delegates to __str__ method for contained objects
        """
        strucStr =  " angle  %s - %s - %s type %s "%(self.ptype1,self.ptype2,self.ptype3,self.type)
        
        if( self.type ==  "harmonic" ):
            strucStr += "\n  harmonic theta_0 = %f K = %f " %(self.theta0 ,self.kb )
        
        return strucStr


    def setharmonic(self, theta0, kb):
        """
        set Harmonic parameters

        E = kb( theta - theta_0 )^2 

        Args:
            theta0 (float) angle             deg  
            kb     (float) force constant    kcal/mol
        """

        if isinstance(theta0, float):
            self.theta0 = theta0
        else:
            print "1st arg should be float"
            raise TypeError

        if isinstance(kb, float):
            self.kb = kb
        else:
            print "2nd arg should be float"
            raise TypeError


class dihtype:
    """
    Set of Dihedral angle parameters 
    """

    def __init__(self, ptype1="blank", ptype2="blank", ptype3="blank", ptype4="blank" , type="harmonic" ):
        """
        Constructor for a angle parameter.
        
        Args:
             ptype1  (str)   Atom type 
             ptype2  (str)   Atom type 
             ptype3  (str)   Atom type 
             ptype4  (str)   Atom type 
             type    (str)   Bond type 
        """

        if isinstance(ptype1, str):
            self.ptype1 = ptype1
        else:
            print "1st arg should be str"
            raise TypeError

        if isinstance(ptype2, str):
            self.ptype2 = ptype2
        else:
            print "2nd arg should be str"
            raise TypeError

        if isinstance(ptype3, str):
            self.ptype3 = ptype3
        else:
            print "3rd arg should be str"
            raise TypeError

        if isinstance(ptype4, str):
            self.ptype4 = ptype4
        else:
            print "4th arg should be str"
            raise TypeError

        if isinstance(type, str):
            self.type = type
        else:
            print "5th arg should be str"
            raise TypeError

        # Set default values for parameters
        self.d = 1.0
        self.mult = 0.0
        self.kb = 0.0
        self.theat_s = 0.0

        # Fourier opls coefficients 
        self.k1 = 0.0 
        self.k2 = 0.0 
        self.k3 = 0.0 
        self.k4 = 0.0 

        # Ryckaert-Bellemans function coefficients 
        self.C0 = 0.0 
        self.C1 = 0.0 
        self.C2 = 0.0 
        self.C3 = 0.0 
        self.C4 = 0.0 
        self.C5 = 0.0 

    def __del__(self):
        """
        Destructor, clears object memory
        """

        
    def __str__(self):
        """
        'Magic' method for printng contents
        Delegates to __str__ method for contained objects
        """
        strucStr =  " dihedral  %s - %s - %s - %s type %s "%(self.ptype1,self.ptype2,self.ptype3,self.ptype4,self.type)
        
        if( self.type ==  "harmonic" ):
            strucStr += "\n  harmonic d = %f n = %f K = %f " %(self.d,self.n ,self.kb )
        
        return strucStr


    def setharmonic(self, d, mult, kb,theat_s):
        """
        set Harmonic parameters

        E = kb[ 1 - d cos( mult theta - theat_s ) ]

        Args:
            d     (float) 
            mult     (float) 
            kb     (float) force constant    kcal/mol
        """

        if isinstance(d, float):
            self.d = d
        else:
            print "1st arg should be float"
            raise TypeError

        if isinstance(kb, float):
            self.kb = kb
        else:
            print "2nd arg should be float"
            raise TypeError

        if isinstance(mult, float):
            self.mult = mult
        else:
            print "3rd arg should be float"
            raise TypeError


    def setopls(self,k1,k2,k3,k4):
        """
        set opls parameters

        E = 1/2 k1[1+cos(theta)]+1/2 k2[1-cos(2 theta)]+1/2 k3[1+cos(3 theta)]+1/2 k4[1-cos(4 theta)]

        Args:
            k1     (float) force constant    kcal/mol
            k2     (float) force constant    kcal/mol
            k3     (float) force constant    kcal/mol
            k4     (float) force constant    kcal/mol
        """

        if isinstance(k1, float):
            self.k1 = k1
        else:
            print "1st arg should be float"
            raise TypeError

        if isinstance(k2, float):
            self.k2 = k2
        else:
            print "2nd arg should be float"
            raise TypeError

        if isinstance(k3, float):
            self.k3 = k3
        else:
            print "3rd arg should be float"
            raise TypeError

        if isinstance(k4, float):
            self.k4 = k4
        else:
            print "4th arg should be float"
            raise TypeError


    def setrb(self,C0,C1,C2,C3,C4,C5):
        """
        set Ryckaert-Bellemans parameters

        V_{rb}(theta) = \sum_n=0^5 C_n [ cos(theata - 180 )^n ]

        Args:
            C0     (float) force constant    kcal/mol
            C1     (float) force constant    kcal/mol
            C2     (float) force constant    kcal/mol
            C3     (float) force constant    kcal/mol
            C4     (float) force constant    kcal/mol
            C5     (float) force constant    kcal/mol
        """

        if isinstance(C0, float):
            self.C0 = C0
        else:
            print "1st arg should be float"
            raise TypeError

        if isinstance(C1, float):
            self.C1 = C1
        else:
            print "2nd arg should be float"
            raise TypeError

        if isinstance(C2, float):
            self.C2 = C2
        else:
            print "3rd arg should be float"
            raise TypeError

        if isinstance(C3, float):
            self.C3 = C3
        else:
            print "4th arg should be float"
            raise TypeError

        if isinstance(C4, float):
            self.C4 = C4
        else:
            print "5th arg should be float"
            raise TypeError

        if isinstance(C5, float):
            self.C5 = C5
        else:
            print "6th arg should be float"
            raise TypeError


class LJtypesContainer:
    """
    Container for ljtypes 
    """

    def __init__(self ):
        """
        Constructor for LJtypesContainer.
        """
        self.ljtypes=dict()                              # Creates empty dict struc

        self.maxgid = 0            # default=0 if idList empty
        
        #if len(idList) == 0:         # If list not set in constructor arg
        #    self.maxgid=0            # default=0 if idList empty
        #else:                        #
        #    self.maxgid=max(idList)  # take max in list for maxgid

    def __del__(self):
        """
        Destructor, clears object memory
        """
    
    def put(self, ljtyp ):
        """
        Append ljtype object to this container.

        Args:
            ljtyp (ljtype) correctly initialized ljtype object

        """
        if isinstance(ljtyp, ljtype):
            self.maxgid += 1
            self.ljtypes[self.maxgid] = copy.deepcopy(ljtyp)
        else:
            print "Attempting to add non-ljtype type to container"
            raise TypeError

    
    def __iter__(self):
        """
        'Magic' method implementing (for x in 'this')....
        """
        # return iter(self.particles)
        return self.ljtypes.iteritems()
    


class BondtypesContainer:
    """
    Container for bondtype
    """

    def __init__(self ):
        """
        Constructor for BondtypesContainer.
        """
        self.bondtypes=dict()                              # Creates empty dict struc

        self.maxgid = 0            # default=0 if idList empty
        
        #if len(idList) == 0:         # If list not set in constructor arg
        #    self.maxgid=0            # default=0 if idList empty
        #else:                        #
        #    self.maxgid=max(idList)  # take max in list for maxgid

    def __del__(self):
        """
        Destructor, clears object memory
        """
    
    def put(self, btyp ):
        """
        Append bondtype object to this container.

        Args:
            btyp (bondtype) correctly initialized bondtype object

        """
        if isinstance(btyp, bondtype):
            self.maxgid += 1
            self.bondtypes[self.maxgid] = copy.deepcopy(btyp)
        else:
            print "Attempting to add non-bondtypes type to container"
            raise TypeError

    
    def __iter__(self):
        """
        'Magic' method implementing (for x in 'this')....
        """
        # return iter(self.particles)
        return self.bondtypes.iteritems()
    


class AngletypesContainer:
    """
    Container for angletype 
    """

    def __init__(self ):
        """
        Constructor for AngletypesContainer.
        """
        self.angletypes=dict()                              # Creates empty dict struc

        self.maxgid = 0            # default=0 if idList empty
        
        #if len(idList) == 0:         # If list not set in constructor arg
        #    self.maxgid=0            # default=0 if idList empty
        #else:                        #
        #    self.maxgid=max(idList)  # take max in list for maxgid

    def __del__(self):
        """
        Destructor, clears object memory
        """
    
    def put(self, atyp ):
        """
        Append angletype object to this container.

        Args:
            atyp (angletype) correctly initialized angletype object

        """
        if isinstance(atyp, angletype):
            self.maxgid += 1
            self.angletypes[self.maxgid] = copy.deepcopy(atyp)
        else:
            print "Attempting to add non-angletypes type to container"
            raise TypeError

    
    def __iter__(self):
        """
        'Magic' method implementing (for x in 'this')....
        """
        # return iter(self.particles)
        return self.angletypes.iteritems()
    

class DihtypesContainer:
    """
    Container for dihtype 
    """

    def __init__(self ):
        """
        Constructor for DihtypesContainer.
        """
        self.dihtypes=dict()                              # Creates empty dict struc

        self.maxgid = 0            # default=0 if idList empty
        
        #if len(idList) == 0:         # If list not set in constructor arg
        #    self.maxgid=0            # default=0 if idList empty
        #else:                        #
        #    self.maxgid=max(idList)  # take max in list for maxgid

    def __del__(self):
        """
        Destructor, clears object memory
        """
    
    def put(self, dtyp ):
        """
        Append dihtype object to this container.

        Args:
            dtyp (dihtype) correctly initialized dihtype object

        """
        if isinstance(dtyp, dihtype):
            self.maxgid += 1
            self.dihtypes[self.maxgid] = copy.deepcopy(dtyp)
        else:
            print "Attempting to add non-dihtype type to container"
            raise TypeError

    
    def __iter__(self):
        """
        'Magic' method implementing (for x in 'this')....
        """
        return self.dihtypes.iteritems()
    

class ParameterContainer:
    """
    Container for force-field parameters
    """

    def __init__(self, ljtypC = LJtypesContainer(), btypC = BondtypesContainer(),atypC=AngletypesContainer(),dtypC=DihtypesContainer() ):
        """
        Constructor for ParameterContainer.

        Args:
             ljtypC  (LJtypesContainer)
             btypC  (BondtypesContainer)
             atypC  (AngletypesContainer)
             dtypC  (DihtypesContainer)
        """
        
        if isinstance(ljtypC, LJtypesContainer):
            self.ljtypC = copy.deepcopy(ljtypC)
        else:
            print "1st arg should be a LJtypesContainer object"
            raise TypeError
        
        if isinstance(btypC, BondtypesContainer):
            self.btypC = copy.deepcopy(btypC)
        else:
            print "2nd arg should be a BondtypesContainer object"
            raise TypeError
        
        if isinstance(atypC, AngletypesContainer):
            self.atypC = copy.deepcopy(atypC)
        else:
            print "3rd arg should be a AngletypesContainer object"
            raise TypeError
        
        if isinstance(dtypC, DihtypesContainer):
            self.dtypC = copy.deepcopy(dtypC)
        else:
            print "4th arg should be a DihtypesContainer object"
            raise TypeError

        
    def __del__(self):
        """
        Destructor, clears object memory
        """

