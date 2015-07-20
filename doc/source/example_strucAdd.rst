
Test for combining StructureContainer-s
=======================================

:download:`Download ipython version <example_strucAdd.ipynb>`

This test shows how to add StructureContainer objects together and how indexing changes within each object
----------------------------------------------------------------------------------------------------------

Load Particle, Bond and StructureContainer modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    #!/usr/bin/env python
    import os, sys, math, random, time
    
    from particles import Particle
    from particles import ParticleContainer
    
    from bonds import Bond
    from bonds import BondContainer
    
    from structureContainer import StructureContainer

Create set of test particle and bond objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    p1 = Particle( [0.2, 1.3,  33.0], "Si", 2.0, 1.23)
    p2 = Particle( [5.0, 2.3, -22.1], "C",  1.0, 2.34)
    p3 = Particle( [5.0, 2.3, -20.1], "C",  1.0, 2.34)
    
    b1 = Bond( 1, 2, 1.233, "hooke")
    b2 = Bond( 2, 3, 0.500, "hooke")

Push particle and bond objects into their respective containers and then use these to build the 'first' StructureContainer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    atoms1   = ParticleContainer()
    atoms1.put(p1)
    atoms1.put(p2)
    atoms1.put(p3)
    
    bonds1   = BondContainer()
    bonds1.put(b1)
    bonds1.put(b2)
    
    polymer1 = StructureContainer(atoms1, bonds1)  # Complete structure 1 completely

.. parsed-literal::

    Cleaning structureContainer


Create Particle, Bond objects. Then create their containers and final the StructureContainer for the 'second' structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    p1other = Particle( [0.0, 2.3, -20.1], "C",  1.0, 2.34)
    p2other = Particle( [50.0, 0.3, -0.1], "Ar", 2.0, 2.34)
    
    b1other = Bond( 1, 2, 1.233, "hooke")    # Correct ptclIDs for second structure
    
    atoms2 = ParticleContainer()
    atoms2.put(p1other)
    atoms2.put(p2other)
    
    bonds2   = BondContainer()
    bonds2.put(b1other)
    
    polymer2 = StructureContainer(atoms2, bonds2)  # Complete structure 2
    print "Number of particles in polymer2 = ", polymer2.getPtclNum()

.. parsed-literal::

    Cleaning structureContainer
    Number of particles in polymer2 =  2


Clean all auxillary data objects, then one is left with only the StructureContainer-s '1' and '2'
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    del p1, p2, p3, p1other, p2other, b1, b2, b1other, atoms1, atoms2, bonds1, bonds2
    print "\n Cleaning memory for initial objects \n" 

.. parsed-literal::

    
     Cleaning memory for initial objects 
    


Structure containers initial state, before adding
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    print "polymer1 = ", polymer1
    print "polymer2 = ", polymer2

.. parsed-literal::

    polymer1 =  
    ---------------------------------------------------------------------
        Structure properties 
    ---------------------------------------------------------------------
          Box lengths: 
            Lx (A) = [0.0, 1.0]
            Ly (A) = [0.0, 1.0]
            Lz (A) = [0.0, 1.0]
          Volume 1000000.000000  A^3 
          Mass 5.910000  AMU 
          Density 0.000010 g/cm^3 
          Lattice vectors 
            v_i (A)  ( 100.000000 , 0.000000 , 0.000000 ) 
            v_j (A)  ( 0.000000 , 100.000000 , 0.000000 ) 
            v_k (A)  ( 0.000000 , 0.000000 , 100.000000 ) 
    
          Particles 3 
          Bonds  2 
          Angles 0 
          Dihedrals 0 
          Impropers 0 
    
     Contains particle objects: 
     1 :  Si 0.200000 1.300000 33.000000 2.000000 1.230000    
     2 :  C 5.000000 2.300000 -22.100000 1.000000 2.340000    
     3 :  C 5.000000 2.300000 -20.100000 1.000000 2.340000    
    
     Contains bond objects: 
     1 :  1 - 2    hooke  
     2 :  2 - 3    hooke  
    
     Contains angle objects: 
    
     Contains dihedral objects: 
    
    polymer2 =  
    ---------------------------------------------------------------------
        Structure properties 
    ---------------------------------------------------------------------
          Box lengths: 
            Lx (A) = [0.0, 1.0]
            Ly (A) = [0.0, 1.0]
            Lz (A) = [0.0, 1.0]
          Volume 1000000.000000  A^3 
          Mass 4.680000  AMU 
          Density 0.000008 g/cm^3 
          Lattice vectors 
            v_i (A)  ( 100.000000 , 0.000000 , 0.000000 ) 
            v_j (A)  ( 0.000000 , 100.000000 , 0.000000 ) 
            v_k (A)  ( 0.000000 , 0.000000 , 100.000000 ) 
    
          Particles 2 
          Bonds  1 
          Angles 0 
          Dihedrals 0 
          Impropers 0 
    
     Contains particle objects: 
     1 :  C 0.000000 2.300000 -20.100000 1.000000 2.340000    
     2 :  Ar 50.000000 0.300000 -0.100000 2.000000 2.340000    
    
     Contains bond objects: 
     1 :  1 - 2    hooke  
    
     Contains angle objects: 
    
     Contains dihedral objects: 
    


Use the '+=' magic method to add the contents of polymer2 into the StructureContainer polymer1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    polymer1 += polymer2
    print "polymer1 = ", polymer1
    print "Number of particles in polymer1 after add = ", polymer1.getPtclNum()

.. parsed-literal::

    polymer1 =  
    ---------------------------------------------------------------------
        Structure properties 
    ---------------------------------------------------------------------
          Box lengths: 
            Lx (A) = [0.0, 1.0]
            Ly (A) = [0.0, 1.0]
            Lz (A) = [0.0, 1.0]
          Volume 1000000.000000  A^3 
          Mass 10.590000  AMU 
          Density 0.000018 g/cm^3 
          Lattice vectors 
            v_i (A)  ( 100.000000 , 0.000000 , 0.000000 ) 
            v_j (A)  ( 0.000000 , 100.000000 , 0.000000 ) 
            v_k (A)  ( 0.000000 , 0.000000 , 100.000000 ) 
    
          Particles 5 
          Bonds  3 
          Angles 0 
          Dihedrals 0 
          Impropers 0 
    
     Contains particle objects: 
     1 :  Si 0.200000 1.300000 33.000000 2.000000 1.230000    
     2 :  C 5.000000 2.300000 -22.100000 1.000000 2.340000    
     3 :  C 5.000000 2.300000 -20.100000 1.000000 2.340000    
     4 :  C 0.000000 2.300000 -20.100000 1.000000 2.340000    
     5 :  Ar 50.000000 0.300000 -0.100000 2.000000 2.340000    
    
     Contains bond objects: 
     1 :  1 - 2    hooke  
     2 :  2 - 3    hooke  
     3 :  4 - 5    hooke  
    
     Contains angle objects: 
    
     Contains dihedral objects: 
    
    Number of particles in polymer1 after add =  5


-------------------- Results (check above) --------------------
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1---b1---2---b2---3 + 1---b1----2 should go to 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1---b1---2---b2---3 4---b3----5 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After adding, polymer2 should be unchanged
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    print "polymer2 = ", polymer2

.. parsed-literal::

    polymer2 =  
    ---------------------------------------------------------------------
        Structure properties 
    ---------------------------------------------------------------------
          Box lengths: 
            Lx (A) = [0.0, 1.0]
            Ly (A) = [0.0, 1.0]
            Lz (A) = [0.0, 1.0]
          Volume 1000000.000000  A^3 
          Mass 4.680000  AMU 
          Density 0.000008 g/cm^3 
          Lattice vectors 
            v_i (A)  ( 100.000000 , 0.000000 , 0.000000 ) 
            v_j (A)  ( 0.000000 , 100.000000 , 0.000000 ) 
            v_k (A)  ( 0.000000 , 0.000000 , 100.000000 ) 
    
          Particles 2 
          Bonds  1 
          Angles 0 
          Dihedrals 0 
          Impropers 0 
    
     Contains particle objects: 
     1 :  C 0.000000 2.300000 -20.100000 1.000000 2.340000    
     2 :  Ar 50.000000 0.300000 -0.100000 2.000000 2.340000    
    
     Contains bond objects: 
     1 :  1 - 2    hooke  
    
     Contains angle objects: 
    
     Contains dihedral objects: 
    

