"""
Class data structures force-field parameters 
"""
import copy, sys

sperator_line = "---------------------------------------------------------------------\n"

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
        self.charge = 0.0 
        self.pid = 0
        self.ptype = "A"  # atomic
        self.lmpindx = 0

    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.ptype1
        del self.epsilon
        del self.sigma
        del self.mass 
        del self.charge
        del self.pid
        del self.ptype 
        
    def __str__(self):
        """
        'Magic' method for printng contents
        Delegates to __str__ method for contained objects
        """
        strucStr =  " LJ  %s  epsilon %f sigma %f "%(self.ptype1, self.epsilon, self.sigma )
                
        return strucStr

    def setpid(self,pid):
        """
        set particle ID

        Args:
            pid (int)  element number for atomic types 
        """

        if isinstance(pid, int):
            self.pid = pid
        else:
            print "1st arg should be int"
            raise TypeError

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
        

    def get_epsilon(self,en_units="kcal"):
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

        # Lammps and gromacs index
        self.lmpindx = 0 
        self.g_indx = 0 


    def set_lmpindx(self,lmpindx):
        """
        Set bond type index for lammps
        """
        self.lmpindx = lmpindx
        
        
    def get_lmpindx(self):
        """
        Return bond type index for lammps
        """
        return self.lmpindx


    def set_g_indx(self,g_indx):
        """
        Set bond type index for gromacs 
        """
        self.g_indx = g_indx
        
        
    def get_g_indx(self):
        """
        Return bond type index for gromacs
        """
        return self.g_indx


    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.ptype1 
        del self.ptype2
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
        strucStr =  " bond  %s - %s type %s "%(self.ptype1,self.ptype2,self.type)
        
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

        # Lammps and gromacs index
        self.lmpindx = 0 
        self.g_indx = 0 


    def set_lmpindx(self,lmpindx):
        """
        Set bond type index for lammps
        """
        self.lmpindx = lmpindx
        
        
    def get_lmpindx(self):
        """
        Return bond type index for lammps
        """
        return self.lmpindx


    def set_g_indx(self,g_indx):
        """
        Set bond type index for gromacs 
        """
        self.g_indx = g_indx
        
        
    def get_g_indx(self):
        """
        Return bond type index for gromacs
        """
        return self.g_indx

    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.ptype1
        del self.ptype2 
        del self.ptype3
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
        strucStr =  " angle  %s - %s - %s type %s "%(self.ptype1,self.ptype2,self.ptype3,self.type)
        
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

class imptype:
    """
    Set of improper dihedral angle parameters
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
        self.e0 = 0.0
        self.ke = 1.0
        self.pn = 0.0    # For periodicimproper

        # Lammps and gromacs index
        self.lmpindx = 0 
        self.g_indx = 0 


    def set_lmpindx(self,lmpindx):
        """
        Set bond type index for lammps
        """
        self.lmpindx = lmpindx
        
        
    def get_lmpindx(self):
        """
        Return bond type index for lammps
        """
        return self.lmpindx


    def set_g_indx(self,g_indx):
        """
        Set bond type index for gromacs 
        """
        self.g_indx = g_indx
        
        
    def get_g_indx(self):
        """
        Return bond type index for gromacs
        """
        return self.g_indx

    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.ptype1
        del self.ptype2 
        del self.ptype3
        del self.ptype4
        del self.ke
        del self.e0
        del self.lmpindx
        del self.g_indx 

        
    def __str__(self):
        """
        'Magic' method for printng contents
        Delegates to __str__ method for contained objects
        """
        strucStr =  " dihedral  %s - %s - %s - %s type %s "%(self.ptype1,self.ptype2,self.ptype3,self.ptype4,self.type)
        
        if( self.type ==  "harmonic" ):
            strucStr += "\n  imp e0 = %f ke = %f lammps index %d  gromcas index %d " %(self.e0,self.ke,self.lmpindx ,self.g_indx )

        return strucStr



    def setimp(self, e0, ke):
        """
        set Harmonic parameters

        E = ? TRAVIS

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

    def getimp(self):
        """
        Get Improper parameters

        Return:
            e0     (float) 
            kb     (float) force constant    kcal/mol
        """
        return self.e0, self.ke

    def set_pn(self, pn):
        """
        set Periodic Harmonic parameter pn 

        Args:
            pn     (int) 
        """

        if isinstance(pn, int):
            self.pn = pn
        else:
            print "1st arg should be int"
            raise TypeError


    def get_pn(self):
        """
        get Periodic Harmonic parameter pn 

        Return:
            pn     (int) 
        """
        return self.np
                    
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


        self.e0 = 0.0
        self.ke = 1.0 

        # Lammps and gromacs index
        self.lmpindx = 0 
        self.g_indx = 0 


    def set_lmpindx(self,lmpindx):
        """
        Set bond type index for lammps
        """

        if isinstance(lmpindx, int ):
            self.lmpindx = lmpindx
        else:
            print "1st arg should be int"
            raise TypeError

        #if( lmpindx <= 0 ):
        #    print " lmpindx less than 0 "
        #    raise TypeError
        
        
    def get_lmpindx(self):
        """
        Return bond type index for lammps
        """
        return self.lmpindx


    def set_g_indx(self,g_indx):
        """
        Set bond type index for gromacs 
        """
        self.g_indx = g_indx
        
        
    def get_g_indx(self):
        """
        Return bond type index for gromacs
        """
        return self.g_indx

    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.ptype1
        del self.ptype2 
        del self.ptype3
        del self.ptype4
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
        strucStr =  " dihedral  %s - %s - %s - %s type %s "%(self.ptype1,self.ptype2,self.ptype3,self.ptype4,self.type)
        
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


    def setimp(self, e0, ke):
        """
        set Harmonic parameters

        E = ? TRAVIS

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

    def getimp(self):
        """
        Get Improper parameters

        Return:
            e0     (float) 
            kb     (float) force constant    kcal/mol
        """
        return self.e0, self.ke
    
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

    def normforceconstants(self,dihen_norm):
        """
        Divide all force constants by a value
        """
        
        # Harmonic
        self.kb = self.kb/dihen_norm

        # Fourier opls coefficients 
        self.k1 = self.k1/dihen_norm
        self.k2 = self.k2/dihen_norm
        self.k3 = self.k3/dihen_norm
        self.k4 = self.k4/dihen_norm

        # Ryckaert-Bellemans function coefficients 
        self.C0 = self.C0/dihen_norm
        self.C1 = self.C1/dihen_norm
        self.C2 = self.C2/dihen_norm
        self.C3 = self.C3/dihen_norm
        self.C4 = self.C4/dihen_norm
        self.C5 = self.C5/dihen_norm

        
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

    def get_theat_s(self):
        """
        Return dih theat_s
        """
        return self.theat_s

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
        del self.ljtypes
        del self.maxgid 

        
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
    
    def __iadd__(self, ljtypC_b ):
        """
        'Magic' method to implement the '+=' operator

        Only add parameters for types not already in parameter container 
        """

        for indx_i, Obj_i in ljtypC_b:
            add_param = True
            for indx_j, Obj_j in self:
                if( Obj_i.ptype1 == Obj_j.ptype1 ):
                    add_param = False
            if( add_param ):
                self.put(Obj_i)
                
                
        return self

    def __len__(self):
        """
        'Magic' method for returning size of container
        """
        return len(self.ljtypes)

    def __getitem__(self, gid):
        """
        'Magic' method implementing obj[] operator.
        Operations on returned elements change container
        """
        return self.ljtypes[gid]

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
        del self.bondtypes
        del self.maxgid
        
    def __iadd__(self, btypC_b ):
        """
        'Magic' method to implement the '+=' operator         

        Only add parameters for types not already in parameter container 
        """

        for indx_i, Obj_i in btypC_b:
            add_param = True
            for indx_j, Obj_j in self:
                if( Obj_i.ptype1 == Obj_j.ptype1 and  Obj_i.ptype2 == Obj_j.ptype2 ):
                    add_param = False
            if( add_param ):
                self.put(Obj_i)
                
        return self

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


    def __getitem__(self, gid):
        """
        'Magic' method implementing obj[] operator.
        Operations on returned elements change container
        """
        return self.bondtypes[gid]

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
        del self.angletypes
        del self.maxgid 

    def __iadd__(self, atypC_b ):
        """
        'Magic' method to implement the '+=' operator
        
        Only add parameters for types not already in parameter container
        
        """

        for indx_i, Obj_i in atypC_b:
            add_param = True
            for indx_j, Obj_j in self:
                if( Obj_i.ptype1 == Obj_j.ptype1 and  Obj_i.ptype2 == Obj_j.ptype2 and  Obj_i.ptype3 == Obj_j.ptype3 ):
                    add_param = False
            if( add_param ):
                self.put(Obj_i)
                
        return self
    
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

    def __getitem__(self, gid):
        """
        'Magic' method implementing obj[] operator.
        Operations on returned elements change container
        """
        return self.angletypes[gid]

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
        del self.dihtypes
        del self.maxgid 

    def __iadd__(self, dtypC_b ):
        """
        'Magic' method to implement the '+=' operator         
        Only add parameters for types not already in parameter container
        
        """

        for indx_i, Obj_i in dtypC_b:
            add_param = True
            for indx_j, Obj_j in self:
                if( Obj_i.ptype1 == Obj_j.ptype1 and  Obj_i.ptype2 == Obj_j.ptype2 and  Obj_i.ptype3 == Obj_j.ptype3 and  Obj_i.ptype4 == Obj_j.ptype4 ):
                    add_param = False
            if( add_param ):
                self.put(Obj_i)
                
        return self
    
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

    def __getitem__(self, gid):
        """
        'Magic' method implementing obj[] operator.
        Operations on returned elements change container
        """
        return self.dihtypes[gid]


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



class ImptypesContainer:
    """
    Container for imptype 
    """

    def __init__(self ):
        """
        Constructor for ImptypesContainer.
        """
        self.imptypes=dict()                              # Creates empty dict struc

        self.maxgid = 0            # default=0 if idList empty
        
        #if len(idList) == 0:         # If list not set in constructor arg
        #    self.maxgid=0            # default=0 if idList empty
        #else:                        #
        #    self.maxgid=max(idList)  # take max in list for maxgid

    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.imptypes
        del self.maxgid 

    def __iadd__(self, ImptypC_b ):
        """
        'Magic' method to implement the '+=' operator         
        Only add parameters for types not already in parameter container
        
        """

        for indx_i, Obj_i in ImptypC_b:
            add_param = True
            for indx_j, Obj_j in self:
                if( Obj_i.ptype1 == Obj_j.ptype1 and  Obj_i.ptype2 == Obj_j.ptype2 and  Obj_i.ptype3 == Obj_j.ptype3 and  Obj_i.ptype4 == Obj_j.ptype4 ):
                    add_param = False
            if( add_param ):
                self.put(Obj_i)
                
        return self
    
    def put(self, Imptyp ):
        """
        Append Imptyp object to this container.

        Args:
            Imptyp (imptype) correctly initialized imptype object

        """
        if isinstance(Imptyp, imptype):
            self.maxgid += 1
            self.imptypes[self.maxgid] = copy.deepcopy(Imptyp)
        else:
            print "Attempting to add non-dihtype type to container"
            raise TypeError

    
    def __iter__(self):
        """
        'Magic' method implementing (for x in 'this')....
        """
        return self.imptypes.iteritems()
    

    def __len__(self):
        """
        'Magic' method for returning size of container
        """
        return len(self.imptypes)

    def __getitem__(self, gid):
        """
        'Magic' method implementing obj[] operator.
        Operations on returned elements change container
        """
        return self.imptypes[gid]


    def findtype(self,fftype_k,fftype_i,fftype_j,fftype_l):
        """
        Find angletypes with a given set of particle types k-i-j
        """

        type_found = False
        cnt_check = 0

        for Imp_all, Imptyp_all  in self:
            all_k = Imptyp_all.ptype1 
            all_i = Imptyp_all.ptype2 
            all_j = Imptyp_all.ptype3
            all_l = Imptyp_all.ptype4

            if ( all_k == fftype_k and  all_i == fftype_i and  all_j == fftype_j and all_l == fftype_l   ):
                type_found = True    
            if ( all_l == fftype_k and all_j == fftype_i and   all_i == fftype_j  and all_k == fftype_l   ):
                type_found = True

            
        
            if( type_found ):
                cnt_check += 1
                print " dih cnt_check ",cnt_check
                print Imptyp_all

class ParameterContainer:
    """
    Container for force-field parameters
    """

    def __init__(self, ljtypC = LJtypesContainer(), btypC = BondtypesContainer(),atypC=AngletypesContainer(),dtypC=DihtypesContainer(),imptypC=ImptypesContainer() ):
        """
        Constructor for ParameterContainer.

        Args:
             ljtypC  (LJtypesContainer)
             btypC  (BondtypesContainer)
             atypC  (AngletypesContainer)
             dtypC  (DihtypesContainer)
             imptypC  (DihtypesContainer)
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
        if isinstance(imptypC, ImptypesContainer):
            self.imptypC = copy.deepcopy(imptypC)
        else:
            print "5th arg should be a ImptypesContainer object"
            raise TypeError
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
        del self.ljtypC
        del self.btypC
        del self.atypC
        del self.dtypC
        del self.imptypC
        
    def __str__(self):
        """
        'Magic' method for printng contents
        """

        strucStr =  "\n"
        strucStr += sperator_line
        strucStr += "    Parameters \n"
        strucStr += sperator_line
        strucStr += "      LJ parameters %d \n"%(len(self.ljtypC))
        strucStr += "      Bond parameters %d \n"%(len( self.btypC))
        strucStr += "      Angle parameters %d \n"%(len( self.atypC))
        strucStr += "      Dihedral parameters %d \n"%(len( self.dtypC))
        strucStr += "      Imporper Dihedral parameters %d \n"%(len( self.imptypC))
        return strucStr

    def __iadd__(self, paramC_b ):
        """
        'Magic' method to implement the '+=' operator 
        
        """
        debug = False 

        # Deep copy parameter container
        #if( len(paramC_b.ljtypC) > 0  ):
        if( debug):

            print "len(self.ljtypC) 0 ",len(self.ljtypC)
            print "len(paramC_b.ljtypC)",len(paramC_b.ljtypC)
            print paramC_b.ljtypC

            for indx_i, Obj_i in self.ljtypC:
                print " a lj object ",indx_i, Obj_i
            for indx_i, Obj_i in paramC_b.ljtypC :
                print " b lj object ",indx_i, Obj_i


            print "len(self.ljtypC) 1 ",len(self.ljtypC)


        self.ljtypC += copy.deepcopy( paramC_b.ljtypC )
        self.btypC += copy.deepcopy( paramC_b.btypC )
        self.atypC += copy.deepcopy( paramC_b.atypC )
        self.dtypC += copy.deepcopy( paramC_b.dtypC )
        self.imptypC += copy.deepcopy( paramC_b.imptypC )

        return self

    def set_nbfunc(self, nbfunc):
        '''
        Set nbfunc
        '''
        if isinstance(nbfunc, int):
            self.nbfunc = nbfunc
        else:
            print "1st arg should be int"
            raise TypeError
                
    def get_nbfunc(self):
        '''
        Get nbfunc
        '''
        return self.nbfunc 


    def set_combmixrule(self,combmixrule):
        '''
        Set combmixrule
        '''
        if isinstance(combmixrule, int):
            self.combmixrule = combmixrule
        else:
            print "1st arg should be int"
            raise TypeError
                
    def get_combmixrule(self):
        '''
        Get combmixrule
        '''
        return self.combmixrule 


    def set_genpairs(self,genpairs):
        '''
        Set genpairs
        '''
        if isinstance(genpairs, str):
            self.genpairs = genpairs
        else:
            print "1st arg should be str"
            raise TypeError
                
    def get_genpairs(self):
        '''
        Get genpairs
        '''
        return self.genpairs 


    def set_fudgeLJ(self,fudgeLJ):
        '''
        Set fudgeLJ
        '''
        if isinstance(fudgeLJ, float):
            self.fudgeLJ = fudgeLJ
        else:
            print "1st arg should be float"
            raise TypeError
                
    def get_fudgeLJ(self):
        '''
        Get fudgeLJ
        '''
        return self.fudgeLJ 


    def set_fudgeQQ(self,fudgeQQ):
        '''
        Set fudgeQQ
        '''
        if isinstance(fudgeQQ, float):
            self.fudgeQQ = fudgeQQ
        else:
            print "1st arg should be float"
            raise TypeError
                
    def get_fudgeQQ(self):
        '''
        Get fudgeQQ
        '''
        return self.fudgeQQ 
