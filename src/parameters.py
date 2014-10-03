"""
Class data structures force-field parameters 
"""
import copy, sys

class ljtype:
    """
    Set of LJ particle parameters 
    """

    def __init__(self, ptype1="blank"):
        """
        Constructor for a Lennard-Jones parameter.
        
        Args:
             ptype1  (str)   Force-field particle type 
        """

        if isinstance(ptype1, str):
            self.ptype1 = ptype1
        else:
            print "1st arg should be str"
            raise TypeError


        # Set default values for parameters
        self.epsilon = 0.0
        self.sigma = 0.0 
        self.mass = 0.0 

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

    def setmass(self,mass):
        """
        set mass of LJ particle 


        Args:
            mass (float)  mass  AMU  
        """

        if isinstance(mass, float):
            self.mass = mass
        else:
            print "1st arg should be float"
            raise TypeError

    def setparam(self,epsilon,sigma):
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

    def get_ptype1(self):
        """
        Return ptype1 of LJ particle
        """
        return self.ptype1
    
    def get_mass(self):
        """
        Return mass of LJ particle
        """
        return self.mass
        

    def get_epsilon(self):
        """
        Return epsilon of LJ particle
        """
        return self.epsilon
        

    def get_sigma(self):
        """
        Return sigma of LJ particle
        """
        return self.sigma
        

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
    
    def get_type(self):
        """
        Return bond type
        """
        return self.type

    def get_ptype1(self):
        """
        Return bond ptype1
        """
        return self.ptype1

    def get_ptype2(self):
        """
        Return bond ptype2
        """
        return self.ptype2

    def get_r0(self):
        """
        Return r0 type
        """
        return self.r0

    def get_kb(self):
        """
        Return kb type
        """
        return self.kb

        
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
        set Harmonic angle parameters

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

    def get_type(self):
        """
        Return angle  type
        """
        return self.type

    def get_ptype1(self):
        """
        Return angle  ptype1
        """
        return self.ptype1

    def get_ptype2(self):
        """
        Return angle  ptype2
        """
        return self.ptype2

    def get_ptype3(self):
        """
        Return angle  ptype3
        """
        return self.ptype3

    def get_theta0(self):
        """
        Return theta0 type
        """
        return self.theta0

    def get_kb(self):
        """
        Return kb type
        """
        return self.kb


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
            strucStr += "\n  harmonic d = %f mult = %f K = %f theat_s = %f " %(self.d,self.mult ,self.kb,self.theat_s )
        
        if( self.type ==  "opls" ):
            strucStr += "\n  k1 = %f k2 = %f k3 = %f k4 = %f " %(self.k1,self.k2,self.k3,self.k4 )

        if( self.type ==  "rb" ):
            strucStr += "\n  C0 = %f  C1 = %f C2 = %f C3 = %f C4 = %f  C5 = %f " %(self.C0,self.C1,self.C2,self.C3,self.C4,self.C5 )

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


    def get_rbClist(self):
        """
        Return list of rb constants 0-5 type
        """
        Clist = []
        Clist.append(self.C0)
        Clist.append(self.C1)
        Clist.append(self.C2)
        Clist.append(self.C3)
        Clist.append(self.C4)
        Clist.append(self.C5)

        return Clist

    def get_oplsklist(self):
        """
        Return list of rb constants 0-5 type
        """
        klist = []
        klist.append(self.k1)
        klist.append(self.k2)
        klist.append(self.k3)
        klist.append(self.k4)

        return klist

    def get_type(self):
        """
        Return dih type
        """
        return self.type

    def get_ptype1(self):
        """
        Return dih ptype1
        """
        return self.ptype1

    def get_ptype2(self):
        """
        Return dih ptype2
        """
        return self.ptype2

    def get_ptype3(self):
        """
        Return dih ptype3
        """
        return self.ptype3

    def get_ptype4(self):
        """
        Return dih ptype4
        """
        return self.ptype4

    def get_d(self):
        """
        Return dih d
        """
        return self.d

    def get_kb(self):
        """
        Return dih kb
        """
        return self.kb

    def get_mult(self):
        """
        Return dih mult
        """
        return self.mult


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
    


    def __len__(self):
        """
        'Magic' method for returning size of container
        """
        return len(self.ljtypes)

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
    


    def __len__(self):
        """
        'Magic' method for returning size of container
        """
        return len(self.bondtypes)


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
    

    def __len__(self):
        """
        'Magic' method for returning size of container
        """
        return len(self.angletypes)


    def findtype(self,fftype_k,fftype_i,fftype_j):
        """
        Find angletypes with a given set of particle types k-i-j
        """

        type_found = False
        cnt_check = 0
        for a_all, atypObj_all  in self:
            all_k = atypObj_all.ptype1 
            all_i = atypObj_all.ptype2 
            all_j = atypObj_all.ptype3

            if( fftype_k == all_k  and fftype_i == all_i  and  fftype_j == all_j ):
                type_found = True 
            if( fftype_j == all_k  and  fftype_i == all_i  and  fftype_k == all_j ):
                type_found = True 

            if( type_found ):
                cnt_check += 1
                print atypObj_all


        
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
    

    def __len__(self):
        """
        'Magic' method for returning size of container
        """
        return len(self.dihtypes)



    def findtype(self,fftype_k,fftype_i,fftype_j,fftype_l):
        """
        Find angletypes with a given set of particle types k-i-j
        """

        type_found = False
        cnt_check = 0

        for d_all, dtypObj_all  in self:
            all_k = dtypObj_all.ptype1 
            all_i = dtypObj_all.ptype2 
            all_j = dtypObj_all.ptype3
            all_l = dtypObj_all.ptype4

            if ( all_k == fftype_k and  all_i == fftype_i and  all_j == fftype_j and all_l == fftype_l   ):
                type_found = True    
            if ( all_l == fftype_k and all_j == fftype_i and   all_i == fftype_j  and all_k == fftype_l   ):
                type_found = True

            
        
            if( type_found ):
                cnt_check += 1
                print " dih cnt_check ",cnt_check
                print dtypObj_all



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


