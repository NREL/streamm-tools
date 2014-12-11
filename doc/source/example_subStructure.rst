
Substructure Test for StructureContainer object
===============================================

This test shows how to set up Structure container with Particle and Bond Containers. Shows how IDs changed in StructureContainer propagate to values set in BondContainer for its held particle ID values. Illustrates how a substructure method can return subgroup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    #!/usr/bin/env python
    import copy
    
    # String diagrams for checking code
    diagramBefore="""
    5--'b1'---2
     -        |
      -       'b2'
       'b4'   |
         -    |
           -  3---'b3'---4
    """
    
    diagramAfter1="""
    5--'b1'---2
    """
    
    diagramAfter2="""
    5--'b1'---2
     -        |
      -       'b2'
       'b4'   |
         -    |
           -  3
    """
Load the Particle and Bond classes and their associated container classes. Load the structureContainer class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    import os, sys, math, random, time
    
    from particles import Particle
    from particles import ParticleContainer
    
    from bonds import Bond
    from bonds import BondContainer
    
    from structureContainer import StructureContainer
Create a set of test Particle and Bond objects. The Bond connectivity is set manually to match the diagramBefore string above.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    p1 = Particle( [0.2, 1.3,  33.0], "Si", 2.0, 1.23)
    p2 = Particle( [5.0, 2.3, -22.1], "C",  1.0, 2.34)
    p3 = Particle( [5.0, 2.3, -20.1], "C",  1.0, 2.34)
    p4 = Particle( [0.0, 2.3, -20.1], "C",  1.0, 2.34)
    p5 = Particle( [0.2, 1.3,  33.0], "Si", 2.0, 1.23)
    
    b1 = Bond( 5, 2, 1.233, "hooke")
    b2 = Bond( 2, 3, 0.500, "hooke")
    b3 = Bond( 3, 4, 2.301, "hooke")
    b4 = Bond( 5, 3, 0.828, "hooke")
Pushing the Particle and Bond objects above into the Particle and Bond Container objects.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    atoms1 = ParticleContainer()
    atoms1.put(p1)
    atoms1.put(p2)
    atoms1.put(p3)
    atoms1.put(p4)
    atoms1.put(p5)
    
    # Example of the 'delete' magic method. Removes the particle with index=1
    del atoms1[1]
    
    bonds = BondContainer()
    bonds.put(b1)
    bonds.put(b2)
    bonds.put(b3)
    bonds.put(b4)
    
    # Removing separate particle and bond objects. This is possible because insertion into the Container objects performs deep copies.
    del p1, p2, p3, p4, b1, b2, b3, b4
Initialize StructureContainer with the Particle and Bond Containers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    polymer1 = StructureContainer(atoms1, bonds)
    del atoms1, bonds
StructureContainer objects implement a python magic method that enables a 'print' statement to output the objects' contents. Default box lengths and the associated lattice vectors are included. At the bottom are the particle, bond, angle, dihedral labels and the connectivity associated with these 2, 3, 4 body interactions.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    print polymer1

.. parsed-literal::

    
    ---------------------------------------------------------------------
        Structure properties 
    ---------------------------------------------------------------------
          Box lengths: 
            Lx (A) = [0.0, 1.0]
            Ly (A) = [0.0, 1.0]
            Lz (A) = [0.0, 1.0]
          Volume 1000000.000000  A^3 
          Mass 8.250000  AMU 
          Density 0.000014 g/cm^3 
          Lattice vectors 
            v_i (A)  ( 100.000000 , 0.000000 , 0.000000 ) 
            v_j (A)  ( 0.000000 , 100.000000 , 0.000000 ) 
            v_k (A)  ( 0.000000 , 0.000000 , 100.000000 ) 
    
          Particles 4 
          Bonds  4 
          Angles 0 
          Dihedrals 0 
          Impropers 0 
    
     Contains particle objects: 
     2 :  C 5.000000 2.300000 -22.100000 1.000000 2.340000    
     3 :  C 5.000000 2.300000 -20.100000 1.000000 2.340000    
     4 :  C 0.000000 2.300000 -20.100000 1.000000 2.340000    
     5 :  Si 0.200000 1.300000 33.000000 2.000000 1.230000    
    
     Contains bond objects: 
     1 :  5 - 2    hooke  
     2 :  2 - 3    hooke  
     3 :  3 - 4    hooke  
     4 :  5 - 3    hooke  
    
     Contains angle objects: 
    
     Contains dihedral objects: 
    


The diagram below should match with the Particle labels and Bond lists above. The 'b1', 'b2' ... are the bond labels. NOTE: the particle positions are not true. This illustrates connectivity and labeling only.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    print diagramBefore

.. parsed-literal::

    
    5--'b1'---2
     -        |
      -       'b2'
       'b4'   |
         -    |
           -  3---'b3'---4
    


Below are example tests of the getSubStructure method in the StructureContainer. All ID's are preserved in the new substructure that is returned. NOTE: the polymer1 object is unchanged after the getSubStructure call
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Example 1: Return substructure containing particle ID's --> [5,2]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

    subpolymer = polymer1.getSubStructure([5,2])
    print subpolymer

.. parsed-literal::

    
    ---------------------------------------------------------------------
        Structure properties 
    ---------------------------------------------------------------------
          Box lengths: 
            Lx (A) = [0.0, 1.0]
            Ly (A) = [0.0, 1.0]
            Lz (A) = [0.0, 1.0]
          Volume 1000000.000000  A^3 
          Mass 3.570000  AMU 
          Density 0.000006 g/cm^3 
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
     2 :  C 5.000000 2.300000 -22.100000 1.000000 2.340000    
     5 :  Si 0.200000 1.300000 33.000000 2.000000 1.230000    
    
     Contains bond objects: 
     1 :  5 - 2    hooke  
    
     Contains angle objects: 
    
     Contains dihedral objects: 
    


.. code:: python

    print "Before ", diagramBefore
    print "After ", diagramAfter1

.. parsed-literal::

    Before  
    5--'b1'---2
     -        |
      -       'b2'
       'b4'   |
         -    |
           -  3---'b3'---4
    
    After  
    5--'b1'---2
    


Example 2: Return substructure containing particle ID's --> [2,3,5]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

    subpolymer = polymer1.getSubStructure([2,3,5])
    print subpolymer

.. parsed-literal::

    Cleaning structureContainer
    
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
          Bonds  3 
          Angles 0 
          Dihedrals 0 
          Impropers 0 
    
     Contains particle objects: 
     2 :  C 5.000000 2.300000 -22.100000 1.000000 2.340000    
     3 :  C 5.000000 2.300000 -20.100000 1.000000 2.340000    
     5 :  Si 0.200000 1.300000 33.000000 2.000000 1.230000    
    
     Contains bond objects: 
     1 :  5 - 2    hooke  
     2 :  2 - 3    hooke  
     4 :  5 - 3    hooke  
    
     Contains angle objects: 
    
     Contains dihedral objects: 
    


.. code:: python

    print "Before ", diagramBefore
    print "After ", diagramAfter2

.. parsed-literal::

    Before  
    5--'b1'---2
     -        |
      -       'b2'
       'b4'   |
         -    |
           -  3---'b3'---4
    
    After  
    5--'b1'---2
     -        |
      -       'b2'
       'b4'   |
         -    |
           -  3
    

