"""
Derived class for structure container
"""

try:
    import simulation as sim
except:
    print "Error: structureContainer module not found"
    print "Check setup.sh path script"
    sys.exit(3)


class SimulationLAMMPS1(sim.Simulation):
    """
    Dervied class implementing input/output methods specific to an application
    """

    def __init__(self, name, verbose=False):
        """
        Constructor for derived class. The base class constructor is called
        explicitly
        """

        # Base class constructor is called
        sim.Simulation.__init__(self, name, verbose)

        if verbose:
            print "Simulation derived class 'LAMMPS1' constructor called"

        self.verbose = verbose
        self.pairCoeffDct = dict()
        self.bondCoeffDct = dict()

    def setCoeffs(self, pairCoeffDct=dict(), bondCoeffDct=dict()):
        """
        Set coefficients for input file write. Specific to this derived class
        
        Args:
            pairCoeffDct (dict) dictionary of potential parameters eg...
                {("Si", "epsilon"):2.30, ("Si", "sigma"):1.0, ("C",  "epsilon"):0.50, ("C",  "sigma"): 0.1 }
            bondCoeffDct (dict) ""
        """
        
        if not isinstance(pairCoeffDct, dict):
            print "setCoeffs: pairCoeffDct should be a python dictionary"
            sys.exit(3)

        if not isinstance(bondCoeffDct, dict):
            print "setCoeffs: bondCoeffDct should be a python dictionary"
            sys.exit(3)

        self.pairCoeffDct = pairCoeffDct
        self.bondCoeffDct = bondCoeffDct


    def __del__(self):
        """
        Destructor, clears object memory
        """

        # Base class destructor is called ?? needed
        sim.Simulation.__del__(self)

        if self.verbose:
            print "Cleaning derived simulation object LAMMPS1"

        del self.pairCoeffDct
        del self.bondCoeffDct


    def writeInput(self, inputName):
        """
        Write out a LAMMPS input data file from all available held
        data (particles, bonds, angles, dihedrals)

        Args:
            inputName    (str)  name of LAMMPS input file to write
        """

        n_atoms = len(self.strucC.ptclC)  # Obtaining particle size from container
        n_bonds = len(self.strucC.bondC)  # " "
        n_angles = 0
        n_dihedrals = 0
        n_impropers = 0

        ptclTypeInfo = self.strucC.ptclC.getTypeInfoDict()  # map of "type":[typeIndex, mass, charge]
        bondTypeInfo = self.strucC.bondC.getTypeInfoDict()  # map of "type":typeIndex

        # Returns map of type,parameter tuple and value
        # SWS: particular to this method
        n_atypes = len(ptclTypeInfo)
        n_btypes = len(bondTypeInfo)
        n_angtypes = 0 
        n_dtypes = 0 
        n_imptypes = 0

        xL = self.strucC.boxLengths[0]
        yL = self.strucC.boxLengths[1]
        zL = self.strucC.boxLengths[2]

        # Open file, write header info
        fileObj = open( inputName, 'w' )
        fileObj.write('LAMMPS Data File \n')
        fileObj.write('\n')
        fileObj.write( "%8d  atoms \n" % n_atoms )
        fileObj.write( "%8d  bonds \n" %  n_bonds )
        fileObj.write( "%8d  angles \n" % n_angles )
        fileObj.write( "%8d  dihedrals \n" %  n_dihedrals )
        fileObj.write( "%8d  impropers \n" % n_impropers  )
        fileObj.write('\n')
        fileObj.write( "%8d  atom types \n" % n_atypes  )
        fileObj.write( "%8d  bond types \n" % n_btypes )
        fileObj.write( "%8d  angle types \n" % n_angtypes )
        fileObj.write( "%8d  dihedral types \n" % n_dtypes )
        fileObj.write( "%8d  improper types \n" % n_imptypes )
        fileObj.write('\n')
        fileObj.write( "%16.8f %16.8f   xlo xhi \n" %  (xL[0] , xL[1] ) )
        fileObj.write( "%16.8f %16.8f   ylo yhi \n" %  (yL[0] , yL[1] ) )
        fileObj.write( "%16.8f %16.8f   zlo zhi \n" %  (zL[0] , zL[1] ) )
        fileObj.write('\n')

        massFormatStr = "%5d %16.8f \n"
        fileObj.write('Masses \n')
        fileObj.write('\n')
        for type, info in ptclTypeInfo.iteritems():
            tIndex = info[0]
            mass   = info[1]
            fileObj.write( massFormatStr % ( tIndex, mass ) )
        fileObj.write('\n')
        
        pairCoeffFormatStr = "%5d %12.6f %12.6f  \n"
        fileObj.write('Pair Coeffs \n')
        fileObj.write('\n')
        for typ in ptclTypeInfo.keys(): # list of types eg ["Si", "C", ..]
            info = ptclTypeInfo[typ]    # map of "type":[typeIndex, mass, charge]
            tIndex = info[0]            # type index for 'typ' (eg "Si")
            epsilon = self.pairCoeffDct[(typ, "epsilon")]
            sigma   = self.pairCoeffDct[(typ, "sigma")]
            fileObj.write( pairCoeffFormatStr % (tIndex, epsilon, sigma  ) )
        fileObj.write('\n')

        bondCoeffFormatStr = "%10d %12.6f %12.6f \n"
        if (n_bonds > 0):
            fileObj.write('Bond Coeffs \n')
            fileObj.write('\n')
            for typ in bondTypeInfo.keys(): # list of types eg ["Si", "C", ..]
                tIndex  = bondTypeInfo[typ]    # map of "type":[typeIndex, mass, charge]
                kenergy = self.bondCoeffDct[(typ, "Kenergy")]
                r0      = self.bondCoeffDct[(typ, "r0")]
                fileObj.write( bondCoeffFormatStr % (tIndex, kenergy, r0) )
            fileObj.write('\n')

        ptclFormatStr = "%5d %5d %5d %12.8f %12.6f %12.6f %12.6f \n"
        fileObj.write('Atoms \n')
        fileObj.write('\n')
        for pid, ptclObj in self.strucC.ptclC:
            pos = ptclObj.position
            mol = ptclObj.tagsDict["molnum"]
            mol = int(mol)
            typ = ptclObj.type
            typeIndex = ptclTypeInfo[typ][0]
            chg = ptclObj.charge
            fileObj.write( ptclFormatStr % (pid, mol, typeIndex, chg, pos[0], pos[1], pos[2] ) )
        fileObj.write('\n')
        
        bondFormatStr = "%9d %8d %9d %9d \n"
        if (n_bonds > 0):
            fileObj.write('Bonds \n')
            fileObj.write('\n')
        for gid, bondObj in self.strucC.bondC:
            pid1 = bondObj.pgid1
            pid2 = bondObj.pgid2
            bondType = 1
            fileObj.write( bondFormatStr % (gid, bondType, pid1, pid2) )
        if (n_bonds > 0):            
            fileObj.write('\n')
        
        # Close LAMMPS file
        fileObj.close()

        #
        # Print type mapping info
        #
        fileObj = open( inputName + ".dat", 'w' )
        fileObj.write("type   typeIndex    mass   charge \n")
        for type, info in ptclTypeInfo.iteritems():
            datFormatStr = "%s %5d %12.6f %12.6f \n"
            tIndex = info[0]
            mass   = info[1]
            charge = info[2]
            fileObj.write(datFormatStr % ( type, tIndex, mass, charge ) )         

        # Close LAMMPS mapping file
        fileObj.close()
