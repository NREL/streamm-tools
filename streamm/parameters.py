"""
parameters module

Class data structures force-field parameters

Most molecular dynamics codes set each particle to certain type (fftype)
usually each interaction is set using these fftypes

"""

import copy 


class LJtype:
    """
    Set of LJ particle parameters 
    """

    def __init__(self, fftype1 ):
        """
        Constructor for a Lennard-Jones parameter.
        
        Args:
             fftype1  (str)   Force-field particle type 
        """

        if isinstance(fftype1, str):
            self.fftype1 = fftype1
        else:
            print "1st arg should be str"
            raise TypeError


        self.properties=dict()
        # Set default values for parameters
        self.epsilon = 0.0
        self.sigma = 0.0 
        self.mass = 0.0 
        self.charge = 0.0 
        self.pkey = 0
        # self.fftype = "A"  # atomic
        self.lmpindx = 0
        self.atomic_symbol = "X"

    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.fftype1
        del self.epsilon
        del self.sigma
        del self.mass 
        del self.charge
        del self.pkey
        del self.lmpindx 
        
    def __str__(self):
        """
        'Magic' method for printng contents
        Delegates to __str__ method for contained objects
        """
        strucStr =  " LJ  %s  epsilon %f sigma %f "%(self.fftype1, self.epsilon, self.sigma )
                
        return strucStr



class Bondtype:
    """
    Set of bond parameters 
    """

    def __init__(self, fftype1='blank', fftype2='blank', type="harmonic" ):
        """
        Constructor for a bond parameter.
        
        Args:
             fftype1  (str)   Atom type 
             fftype2  (str)   Atom type 
             type    (str)   Bond type 
        """

        if isinstance(fftype1, str):
            self.fftype1 = fftype1
        else:
            print "1st arg should be str"
            raise TypeError

        if isinstance(fftype2, str):
            self.fftype2 = fftype2
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

        # Lammps and gromacs index
        self.lmpindx = 0 
        self.g_indx = 0 


    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.fftype1 
        del self.fftype2
        del self.type 
        del self.r0
        del self.kb 
        del self.lmpindx 
        del self.g_indx 
        
    def __str__(self):
        """
        'Magic' method for printng contents
        Delegates to __str__ method for contained objects
        """
        strucStr =  " bond  %s - %s type %s "%(self.fftype1,self.fftype2,self.type)
        
        if( self.type ==  "harmonic" ):
            strucStr += "\n  harmonic r_0 = %f K = %f lammps index %d  gromcas index %d  " %(self.r0 ,self.kb,self.lmpindx ,self.g_indx )
                
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
    
   
class Angletype:
    """
    Set of Angle parameters 
    """

    def __init__(self, fftype1="blank", fftype2="blank", fftype3="blank" , type="harmonic" ):
        """
        Constructor for a angle parameter.
        
        Args:
             fftype1  (str)   Atom type 
             fftype2  (str)   Atom type 
             fftype3  (str)   Atom type 
             type    (str)   Bond type 
        """

        if isinstance(fftype1, str):
            self.fftype1 = fftype1
        else:
            print "1st arg should be str"
            raise TypeError

        if isinstance(fftype2, str):
            self.fftype2 = fftype2
        else:
            print "2nd arg should be str"
            raise TypeError

        if isinstance(fftype3, str):
            self.fftype3 = fftype3
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

        # Lammps and gromacs index
        self.lmpindx = 0 
        self.g_indx = 0 


    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.fftype1
        del self.fftype2 
        del self.fftype3
        del self.type
        del self.theta0
        del self.kb 
        del self.lmpindx
        del self.g_indx 


    def __str__(self):
        """
        'Magic' method for printng contents
        Delegates to __str__ method for contained objects
        """
        strucStr =  " angle  %s - %s - %s type %s "%(self.fftype1,self.fftype2,self.fftype3,self.type)
        
        if( self.type ==  "harmonic" ):
            strucStr += "\n  harmonic theta_0 = %f K = %f lammps index %d  gromcas index %d  " %(self.theta0 ,self.kb,self.lmpindx ,self.g_indx )

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


class Dihtype:
    """
    Set of Dihedral angle parameters 
    """

    def __init__(self, fftype1="blank", fftype2="blank", fftype3="blank", fftype4="blank" , type="multiharmonic" ):
        """
        Constructor for a angle parameter.
        
        Args:
             fftype1  (str)   Atom type 
             fftype2  (str)   Atom type 
             fftype3  (str)   Atom type 
             fftype4  (str)   Atom type 
             type    (str)   Bond type 
        """

        if isinstance(fftype1, str):
            self.fftype1 = fftype1
        else:
            print "1st arg should be str"
            raise TypeError

        if isinstance(fftype2, str):
            self.fftype2 = fftype2
        else:
            print "2nd arg should be str"
            raise TypeError

        if isinstance(fftype3, str):
            self.fftype3 = fftype3
        else:
            print "3rd arg should be str"
            raise TypeError

        if isinstance(fftype4, str):
            self.fftype4 = fftype4
        else:
            print "4th arg should be str"
            raise TypeError

        if isinstance(type, str):
            self.type = type
        else:
            print "5th arg should be str"
            raise TypeError

        # Set default values for parameters
        self.paths = 1 
        self.d = 1.0    # gamma 
        self.mult = 0.0 # n 
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


        self.e0 = 0.0
        self.ke = 1.0 

        # Lammps and gromacs index
        self.lmpindx = 0 
        self.g_indx = 0

    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.fftype1
        del self.fftype2 
        del self.fftype3
        del self.fftype4
        del self.type
        del self.d
        del self.mult
        del self.kb
        del self.theat_s
        del self.k1
        del self.k2
        del self.k3
        del self.k4
        del self.C0
        del self.C1
        del self.C2
        del self.C3
        del self.C4
        del self.C5
        del self.e0
        del self.ke
        del self.lmpindx
        del self.g_indx 

    def __str__(self):
        """
        'Magic' method for printng contents
        Delegates to __str__ method for contained objects
        """
        strucStr =  " dihedral  %s - %s - %s - %s type %s "%(self.fftype1,self.fftype2,self.fftype3,self.fftype4,self.type)
        
        if( self.type ==  "harmonic" ):
            strucStr += "\n  harmonic d = %f mult = %f K = %f theat_s = %f lammps index %d  gromcas index %d " %(self.d,self.mult ,self.kb,self.theat_s,self.lmpindx ,self.g_indx )
        if( self.type ==  "multiharmonic" ):
            strucStr += "\n  harmonic d = %f mult = %f K = %f theat_s = %f lammps index %d  gromcas index %d " %(self.d,self.mult ,self.kb,self.theat_s,self.lmpindx ,self.g_indx )
        if( self.type ==  "opls" ):
            strucStr += "\n  k1 = %f k2 = %f k3 = %f k4 = %f lammps index %d  gromcas index %d " %(self.k1,self.k2,self.k3,self.k4,self.lmpindx ,self.g_indx )
        if( self.type ==  "rb" ):
            strucStr += "\n  C0 = %f  C1 = %f C2 = %f C3 = %f C4 = %f  C5 = %f lammps index %d  gromcas index %d " %(self.C0,self.C1,self.C2,self.C3,self.C4,self.C5 ,self.lmpindx ,self.g_indx)

        return strucStr


    def setharmonic(self, mult, kb,theat_s):
        """
        set MultiHarmonic parameters
        dihedral_style charmm

        E = kb[ 1 - cos( mult theta - theat_s ) ]  gromacs
        E = kb[ 1 - cos( n theta - d ) ]           lammps 

        Args:
            mult     (float) 
            kb     (float) force constant    kcal/mol
            theat_s     (float) angle degrees 
        """

        if isinstance(kb, float):
            self.kb = kb
        else:
            print "1nd arg kb should be float"
            raise TypeError

        if isinstance(mult, float):
            self.mult = mult
        else:
            print "2rd arg mult should be float"
            raise TypeError

        if isinstance(theat_s, float):
            self.theat_s = theat_s
        else:
            print "3th arg theat_s should be float"
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

        # Translate to Ryckaert-Bellemans function
        self.C0 = k2 + 0.5*(k1+k3)
        self.C1 = 0.5*(-1.0*k1+3.0*k3)
        self.C2 = -1.0*k2 + 4.0*k3
        self.C3 = -2.0*k3
        self.C4 = -4.0*k4
        self.C5 = 0.0


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
        # Translate to opls 
        self.k1 = -1.0*( 2.0*C1 + 3.0*C3/2.0)
        self.k2 = -1.0*( C2 + C4)
        self.k3 = -0.5*C3
        self.k4 = -0.25*C4

        
class Imptype:
    """
    Set of improper dihedral angle parameters
    """

    def __init__(self, fftype1="blank", fftype2="blank", fftype3="blank", fftype4="blank" , type="improper" ):
        """
        Constructor for a angle parameter.
        
        Args:
             fftype1  (str)   Atom type 
             fftype2  (str)   Atom type 
             fftype3  (str)   Atom type 
             fftype4  (str)   Atom type 
             type    (str)   Improper type 
        """

        if isinstance(fftype1, str):
            self.fftype1 = fftype1
        else:
            print "1st arg should be str"
            raise TypeError

        if isinstance(fftype2, str):
            self.fftype2 = fftype2
        else:
            print "2nd arg should be str"
            raise TypeError

        if isinstance(fftype3, str):
            self.fftype3 = fftype3
        else:
            print "3rd arg should be str"
            raise TypeError

        if isinstance(fftype4, str):
            self.fftype4 = fftype4
        else:
            print "4th arg should be str"
            raise TypeError

        if isinstance(type, str):
            self.type = type
        else:
            print "5th arg should be str"
            raise TypeError

        # Set default values for parameters
        self.e0 = 0.0
        self.ke = 1.0
        self.pn = 0.0    # For periodicimproper

        # Lammps and gromacs index
        self.lmpindx = 0 
        self.g_indx = 0 

    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.fftype1
        del self.fftype2 
        del self.fftype3
        del self.fftype4
        del self.ke
        del self.e0
        del self.pn
        del self.lmpindx
        del self.g_indx 

        
    def __str__(self):
        """
        'Magic' method for printng contents
        Delegates to __str__ method for contained objects
        """
        strucStr =  " improper  %s - %s - %s - %s type %s "%(self.fftype1,self.fftype2,self.fftype3,self.fftype4,self.type)
        
        if( self.type ==  "improper" ):
            strucStr += "\n  imp e0 = %f ke = %f lammps index %d  gromcas index %d " %(self.e0,self.ke,self.lmpindx ,self.g_indx )

        return strucStr        


    def setimp(self, e0, ke):
        """
        set Harmonic parameters

        E = kb ( e_{lijk} - e0 )^2 

        Args:
            e0     (float) 
            kb     (float) force constant    kcal/mol
        """

        if isinstance(e0, float):
            self.e0 = e0
        else:
            print "1st arg should be float"
            raise TypeError

        if isinstance(ke, float):
            self.ke = ke
        else:
            print "2nd arg should be float"
            raise TypeError



class Container():
    """
    Container for force-field parameters
    """

    def __init__(self, verbose=False):
        """
        Constructor for Parameter.Container

        """
    

        self.ljtypes = dict()                           # Creates empty dict struc
        self.bondtypes = dict()                                   # Creates empty dict struc
        self.angletypes = dict()                                # Creates empty dict struc
        self.dihtypes = dict()                                # Creates empty dict struc
        self.imptypes = dict()                                  # Creates empty dict struc

        # Int count of the length of each dictionary
        #   mostly for internal use 
        self.n_ljtypes = 0    
        self.n_bondtypes = 0    
        self.n_angletypes = 0    
        self.n_dihtypes = 0    
        self.n_imptypes = 0            
        #
        # Set defaults 
        #
        self.nbfunc = 1      # Use 1 (Lennard-Jones) or 2 (Buckingham)
        self.combmixrule = 3   # 1 or 3  - geometric; 2 - arithmetic
        self.genpairs = "yes"  # generates 1-4 parameters that are not present in the pair list 
        self.fudgeLJ = 1.0   # the factor by which to multiply Lennard-Jones 1-4 interactions, default 1
        self.fudgeQQ = 1.0   # the factor by which to multiply electrostatic 1-4 interactions, default 1
        

    def __del__(self):
        """
        Destructor, clears object memory
        """
        # dictionary 
        del self.ljtypes
        del self.bondtypes
        del self.angletypes
        del self.dihtypes
        del self.imptypes
        # 
        del self.nbfunc
        del self.combmixrule
        del self.genpairs
        del self.fudgeLJ
        del self.fudgeQQ
        
    def __str__(self):
        """
        'Magic' method for printng contents
        """
        sperator_line = "---------------------------------------------------------------------\n"

        strucStr =  "\n"
        strucStr += sperator_line
        strucStr += "    Parameters \n"
        strucStr += sperator_line
        strucStr += "      LJ parameters %d \n"%(self.n_ljtypes)
        strucStr += "      Bond parameters %d \n"%(self.n_bondtypes)
        strucStr += "      Angle parameters %d \n"%(self.n_angletypes)
        strucStr += "      Dihedral parameters %d \n"%(self.n_dihtypes)
        strucStr += "      Imporper Dihedral parameters %d \n"%(self.n_imptypes)
        return strucStr

    def __iadd__(self, other ):
        """
        'Magic' method to implement the '+=' operator 
        
        """

        if isinstance(other, Container):
            for ljtkey_i, ljtype_i  in other.ljtypes.iteritems():
                self.add_LJtype(ljtype_i,deepcopy = True )
            for btkey_i,bondtype_i  in other.bondtypes.iteritems():
                self.add_bondtype(bondtype_i,deepcopy = True )
            for atkey_i,angletype_i  in other.angletypes.iteritems():
                self.add_angletype(angletype_i,deepcopy = True )
            for dtkey_i, dihtype_i  in other.dihtypes.iteritems():    
                self.add_dihtype(dihtype_i,deepcopy = True )
            for itkey_i, imptype_i  in other.imptypes.iteritems():    
                self.add_imptype(imptype_i,deepcopy = True )


        else:
            print "1st arg should be a Parameter.Container object"
            raise TypeError
        

        return self


        
    def add_LJtype(self, ljtype_i, deepcopy = True ):
        """
        Add 'Ljtype' object to ljtypes dict in this container and update n_ljtypes accordingly
        """
        if isinstance(ljtype_i, LJtype):
            self.n_ljtypes = len(self.ljtypes)
            if( deepcopy ):
                self.ljtypes[self.n_ljtypes] = copy.deepcopy(ljtype_i) # index 0 -> (N-1)
            else:
                self.ljtypes[self.n_ljtypes] = ljtype_i # index 0 -> (N-1)
                
            self.n_ljtypes = len(self.ljtypes)
        else:
            print "Attempting to add non-Ljtype type to container"
            raise TypeError




    def add_bondtype(self, bondtype_i, deepcopy = True ):
        """
        Add 'Bondtype' object to bondtypes dict in this container and update n_bondtypes accordingly
        """
        if isinstance(bondtype_i, Bondtype):
            self.n_bondtypes = len(self.bondtypes)
            if( deepcopy ):
                self.bondtypes[self.n_bondtypes] = copy.deepcopy(bondtype_i) # index 0 -> (N-1)
            else:
                self.bondtypes[self.n_bondtypes] = bondtype_i # index 0 -> (N-1)
                
            self.n_bondtypes = len(self.bondtypes)
        else:
            print "Attempting to add non-Bondtype type to container"
            raise TypeError


    def add_angletype(self, angletype_i, deepcopy = True ):
        """
        Add 'Angletype' object to angletypes dict in this container and update n_angletypes accordingly
        """
        if isinstance(angletype_i, Angletype):
            self.n_angletypes = len(self.angletypes)
            if( deepcopy ):
                self.angletypes[self.n_angletypes] = copy.deepcopy(angletype_i) # index 0 -> (N-1)
            else:
                self.angletypes[self.n_angletypes] = angletype_i # index 0 -> (N-1)
                
            self.n_angletypes = len(self.angletypes)
        else:
            print "Attempting to add non-Angletype type to container"
            raise TypeError

    def add_dihtype(self, dihtype_i, deepcopy = True ):
        """
        Add 'Dihtype' object to dihtypes dict in this container and update n_dihtypes accordingly
        """
        if isinstance(dihtype_i, Dihtype):
            self.n_dihtypes = len(self.dihtypes)
            if( deepcopy ):
                self.dihtypes[self.n_dihtypes] = copy.deepcopy(dihtype_i) # index 0 -> (N-1)
            else:
                self.dihtypes[self.n_dihtypes] = dihtype_i # index 0 -> (N-1)
                
            self.n_dihtypes = len(self.dihtypes)
        else:
            print "Attempting to add non-Dihtype type to container"
            raise TypeError


    def add_imptype(self, imptype_i, deepcopy = True ):
        """
        Add 'Imptype' object to imptypes dict in this container and update n_imptypes accordingly
        """
        if isinstance(imptype_i, Imptype):
            self.n_imptypes = len(self.imptypes)
            if( deepcopy ):
                self.imptypes[self.n_imptypes] = copy.deepcopy(imptype_i) # index 0 -> (N-1)
            else:
                self.imptypes[self.n_imptypes] = imptype_i # index 0 -> (N-1)
                
            self.n_imptypes = len(self.imptypes)
        else:
            print "Attempting to add non-Imptype type to container"
            raise TypeError
