"""
Class data structures for 2, 3, 4 point groupings of Particle objects
"""
import copy, sys

class Bond:
    """
    Data structure for describing any 2-point associatiaon of Particle-s
    """

    def __init__(self, pgid1=0, pgid2=0, length=0.0, type="blank"):
        """
        Constructor for a general bond. Checks for types in arguments
        and throws a TypeError when appropriate

        Args:
            pgid1   (int)   GlobalID of Particle object in bond
            pgid2   (int)   GlobalID of Particle object in bond
            length (float) Cartesian length of bond
            type   (str)   Charge value in units of [e]
        """
        
        if isinstance(pgid1, int):
            self.pgid1 = pgid1
        else:
            print "1st arg should be int"
            raise TypeError

        if isinstance(pgid2, int):
            self.pgid2 = pgid2
        else:
            print "2nd arg should be int type"
            raise TypeError

        if isinstance(length, float):
            self.length = length
        else:
            print "3rd arg should be float value"
            raise TypeError

        if isinstance(type, str):
            self.type = type
        else:
            print "4th arg should be string value"
            raise TypeError


    def __contains__(self, pgid):
        """
        'Magic' method implementing 'in' keyword.

        Args:
            pgid (int) Particle GID to check against 'held' IDs
        """
        if ( (pgid == self.pgid1) or (pgid == self.pgid2) ):
            return True
        else:
            return False



class BondContainer:
    """
    Main data structure for holding Bond objects. Map of global
    bond ID (integer) to Bond object instances
    """

    def __init__(self):
        """
        Constructor: sets up a dictionary for indexing 'Particle' objects
        """
        self.bonds=dict()
        self.maxgid=0

    def __del__(self):
        """
        Destructor, clears dictionary memory
        """
        del self.bonds


    def __len__(self):
        """
        'Magic' method for returning size of container
        """
        return len(self.bonds)


    def __str__(self):
        """
        'Magic' method for printng contents
        """

        bondStr="\n Contains bond objects: \n"

        for gid in self.bonds:
            bondStr = bondStr + str(gid) + " " + str(self.bonds[gid].__dict__) + "\n"
        return bondStr


    def __setitem__(self, gid, bond):
        """
        'Magic' method implementing obj[]=value operator
        Need
        """
        if gid in self.bonds.keys():
            self.bonds[gid]=copy.deepcopy(bond)
        else:
            print "Cannot add bond object to non-existent ID"
            sys.exit(3) 


    def __getitem__(self, gid):
        """
        'Magic' method implementing obj[] operator
        """
        return self.bonds[gid]


    def __delitem__(self, gid):
        """
        'Magic' method implementing del obj[] operator
        """
        del self.bonds[gid]


    def __iter__(self):
        """
        'Magic' method implementing (for x in 'this')....
        """
        return iter(self.bonds)


    def __contains__(self, gid):
        """
        'Magic' method implementing in keyword (key in obj')....
        """
        return gid in self.bonds


    def __iadd__(self, other):
        """
        'Magic' method to implement the '+=' operator
        
        Note: SWS: not sure what this should mean !?
        """
        return self


    def put(self, bond):
        """
        Append 'Bond' object to this container. Updates globalID for container
        by incrementing the maxgid member

        Args:
            ptcl (Particle) correctly initialized Particle object

        NOTE:
            (1) One can imagine extra conditions on bonds inserted
            (2) This could check for uniqueness of all globalID's and throw error for copies
        """
        
        if isinstance(bond, Bond):
            self.maxgid += 1
            self.bonds[self.maxgid] = copy.deepcopy(bond)
        else:
            print "Attempting to add non-Bond type to container"
            raise TypeError


    def replacePtclIDs(self, findPtclID, newPtclID):
        """
        Replace s that contain globalID of particle

        Args:
            findPtclID (int) globalID of particle to search for
            newPtclID  (int) globalID to replace findPtclID with

        Returns:
        """

        for gid in self.bonds:
            
            bond = self.bonds[gid]
            if findPtclID in bond:

                if bond.pgid1 == findPtclID:
                    bond.pgid1 = newPtclID
                elif bond.pgid2 == findPtclID:
                    bond.pgid2 = newPtclID
                else:
                    print "ptclID not found in bond"
                    sys.exit(3)
