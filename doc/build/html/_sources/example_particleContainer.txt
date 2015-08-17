
Basic Particle class method example
===================================

:download:`Download ipython version <example_particleContainer.ipynb>`

Shows various operators within Particle and ParticleContainer classes. Illustrates memory management structure and access methods.
----------------------------------------------------------------------------------------------------------------------------------

Load Particle modules containing class specifications
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    #!/usr/bin/env python
    
    import os, sys, math, random, time
    
    from particles import Particle
    from particles import ParticleContainer

Create set of test particle data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    p1 = Particle( [0.2, 1.3,  33.0], "Si", 2.0, 1.23)
    p2 = Particle( [5.0, 2.3, -22.1], "C",  1.0, 2.34)
    p3 = Particle( [5.0, 2.3, -20.1], "C",  1.0, 2.34)
    p4 = Particle( [0.0, 2.3, -20.1], "C",  1.0, 2.34)

Push Particle objects into ParticleContainer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    atoms1 = ParticleContainer()
    atoms2 = ParticleContainer()
    
    atoms1.put(p1)
    atoms1.put(p2)
    #
    atoms2.put(p3)
    atoms2.put(p4)
    
    del p1, p2, p3, p4
    print "\n Cleaning memory for initial objects \n" 

.. parsed-literal::

    
     Cleaning memory for initial objects 
    


Assignment operator returns x as a reference to internal data, which can be used to make direct edits on the ParticleContainer data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    print "x = atoms1[1] returns x as an effective 'reference' \n"
    x = atoms1[1]
    print "x = ", x.__dict__, "\n"
    x.position=[1.0, 1.0, 1.0]
    print "after changing with x.position = [1.0, 1.0, 1.0]"
    print "x = ", x.__dict__, "\n"

.. parsed-literal::

    x = atoms1[1] returns x as an effective 'reference' 
    
    x =  {'position': [0.2, 1.3, 33.0], 'charge': 2.0, 'type': 'Si', 'mass': 1.23, 'tagsDict': {'type': 'Si'}} 
    
    after changing with x.position = [1.0, 1.0, 1.0]
    x =  {'position': [1.0, 1.0, 1.0], 'charge': 2.0, 'type': 'Si', 'mass': 1.23, 'tagsDict': {'type': 'Si'}} 
    


This value has been changed by code above
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    print "atoms1 has been changed"
    print atoms1

.. parsed-literal::

    atoms1 has been changed
    
     Contains particle objects: 
     1 :  Si 1.000000 1.000000 1.000000 2.000000 1.230000    
     2 :  C 5.000000 2.300000 -22.100000 1.000000 2.340000    
    


Testing the 'delete' method for ParticleContainer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    print "before, atoms1--> ", atoms1, "\n"
    del atoms1[2]
    print "after 'del atoms1[2]' atoms1 --> ", atoms1, "\n"

.. parsed-literal::

    before, atoms1-->  
     Contains particle objects: 
     1 :  Si 1.000000 1.000000 1.000000 2.000000 1.230000    
     2 :  C 5.000000 2.300000 -22.100000 1.000000 2.340000    
     
    
    after 'del atoms1[2]' atoms1 -->  
     Contains particle objects: 
     1 :  Si 1.000000 1.000000 1.000000 2.000000 1.230000    
     
    


Demonstrates the 'in' operator which can be used to test if an object contains a data member
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    print "Testing 'in' operator (1 in atoms1)"
    if (1 in atoms1):
        print "atoms1 contains gid 1"
    else:
        print "key not found in atoms1"
    
    print " "
        
    print "Testing 'in' operator (5 in atoms1)"
    if (5 in atoms1):
        print "atoms1 contains gid 5"
    else:
        print "key not found in atoms1"

.. parsed-literal::

    Testing 'in' operator (1 in atoms1)
    atoms1 contains gid 1
     
    Testing 'in' operator (5 in atoms1)
    key not found in atoms1


ParticleContainers can be directly combined with the '+=' operator. Note: Particle object indices within a container are always unique. If indices overlap, the maximum index within the container on the left of the '+=' is used to shift all other indices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    print " "
    atoms1 += atoms2
    print "Will print the new atoms1 after adding atoms1 += atoms2"
    print atoms1

.. parsed-literal::

     
    Will print the new atoms1 after adding atoms1 += atoms2
    
     Contains particle objects: 
     1 :  Si 1.000000 1.000000 1.000000 2.000000 1.230000    
     3 :  C 5.000000 2.300000 -20.100000 1.000000 2.340000    
     4 :  C 0.000000 2.300000 -20.100000 1.000000 2.340000    
    

