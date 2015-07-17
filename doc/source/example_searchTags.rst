
Search by Particle tags example
===============================

:download:`Download ipython version <example_searchTags.ipynb>`

This test illustrates the search capability for multiple tags and can be combined with class method that returns iterator over search results
---------------------------------------------------------------------------------------------------------------------------------------------

Load Particle class modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    #!/usr/bin/env python
    import copy
    import os, sys, math, random, time
    
    from particles import Particle
    from particles import ParticleContainer

Create test particle data. Custom tags are also created and set for each Particle object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    p1 = Particle( [0.2, 1.3,  33.0], "Si", 2.0, 1.23)
    tagsD = {"molnum":1,"ringnum":4}
    p1.setTagsDict(tagsD)
    
    p2 = Particle( [5.0, 2.3, -22.1], "C",  1.0, 2.34)
    tagsD = {"molnum":2,"ringnum":4}
    p2.setTagsDict(tagsD)
    
    p3 = Particle( [5.0, 2.3, -20.1], "C",  1.0, 2.34)
    tagsD = {"molnum":1, "ringnum":4}
    p3.setTagsDict(tagsD)
    
    p4 = Particle( [0.0, 2.3, -20.1], "Si",  1.0, 2.34)
    tagsD = {"molnum":2,"ringnum":4}
    p4.setTagsDict(tagsD)
    
    p5 = Particle( [1.0, 2.3, -20.1], "C",  1.0, 5.34)
    tagsD = {"molnum":2,"ringnum":2}
    p5.setTagsDict(tagsD)
    
    p6 = Particle( [8.0, 2.3, -20.1], "Si",  1.0, 8.00)
    tagsD = {"molnum":2,"ringnum":3}
    p6.setTagsDict(tagsD)
    
    p7 = Particle( [8.0, 2.3, -20.1], "O",  1.0, 8.00)
    tagsD = {"molnum":2,"ringnum":3}
    p7.setTagsDict(tagsD)
    
    p8 = Particle( [8.0, 2.3, -20.1], "O",  1.0, 8.00)
    tagsD = {"molnum":2,"ringnum":3}
    p8.setTagsDict(tagsD)
    
    p9 = Particle( [8.0, 2.3, -20.1], "H",  1.0, 8.00)
    tagsD = {"molnum":2,"ringnum":3}
    p9.setTagsDict(tagsD)

Add Particle objects to ParticleContainer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    atoms1 = ParticleContainer()
    atoms1.put(p1)
    atoms1.put(p2)
    atoms1.put(p3)
    atoms1.put(p4)
    atoms1.put(p5)
    atoms1.put(p6)
    atoms1.put(p7)
    atoms1.put(p8)
    atoms1.put(p9)
    del p1, p2, p3, p4, p5, p6, p7, p8, p9
    
    print "atoms1 initially contains all of the following"

.. parsed-literal::

    atoms1 initially contains all of the following


Different search tag combinations Note: default tag like 'Si' can be mixed with user defined tags like 'molnum'
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    searchD = {'molnum':2, 'type':"Si"}
    print "Print only atoms1 with search....", searchD, "\n"
    subList = atoms1.getParticlesWithTags(searchD)
    for id, ptcl in atoms1(subList):
        print "id = ", id, " ptclObj = ",str(ptcl)," molnum ", ptcl.tagsDict["molnum"],"ringnum", ptcl.tagsDict["ringnum"]

.. parsed-literal::

    Print only atoms1 with search.... {'molnum': 2, 'type': 'Si'} 
    
    id =  4  ptclObj =   Si 0.000000 2.300000 -20.100000 1.000000 2.340000    molnum  2 ringnum 4
    id =  6  ptclObj =   Si 8.000000 2.300000 -20.100000 1.000000 8.000000    molnum  2 ringnum 3


.. code:: python

    searchD = {'molnum':2, 'type':"Si",'ringnum':4}
    print "Print only atoms1 with search....", searchD, "\n"
    subList = atoms1.getParticlesWithTags(searchD)
    for id, ptcl in atoms1(subList):
        print "id = ", id, " ptclObj = ",str(ptcl)," molnum ", ptcl.tagsDict["molnum"],"ringnum", ptcl.tagsDict["ringnum"]

.. parsed-literal::

    Print only atoms1 with search.... {'ringnum': 4, 'molnum': 2, 'type': 'Si'} 
    
    id =  4  ptclObj =   Si 0.000000 2.300000 -20.100000 1.000000 2.340000    molnum  2 ringnum 4


Note: multiple values for the same tag can be search (eg 'ringnum) in a list. This example searches for all particles with type = 'Si' and ringnum '3' or '4'
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    searchD = {'type':"Si", 'ringnum':[3, 4]}
    print "Print only atoms1 with search....", searchD, "\n"
    subList = atoms1.getParticlesWithTags(searchD)
    for id, ptcl in atoms1(subList):
        print "id = ", id, " ptclObj = ",str(ptcl)," molnum ", ptcl.tagsDict["molnum"],"ringnum", ptcl.tagsDict["ringnum"]

.. parsed-literal::

    Print only atoms1 with search.... {'ringnum': [3, 4], 'type': 'Si'} 
    
    id =  1  ptclObj =   Si 0.200000 1.300000 33.000000 2.000000 1.230000    molnum  1 ringnum 4
    id =  4  ptclObj =   Si 0.000000 2.300000 -20.100000 1.000000 2.340000    molnum  2 ringnum 4
    id =  6  ptclObj =   Si 8.000000 2.300000 -20.100000 1.000000 8.000000    molnum  2 ringnum 3


.. code:: python

    searchD = {'type':["Si","O"], 'ringnum':[2, 3, 4]}
    print "Print only atoms1 with search....", searchD, "\n"
    subList = atoms1.getParticlesWithTags(searchD)
    for id, ptcl in atoms1(subList):
        print "id = ", id, " ptclObj = ",str(ptcl)," molnum ", ptcl.tagsDict["molnum"],"ringnum", ptcl.tagsDict["ringnum"]

.. parsed-literal::

    Print only atoms1 with search.... {'ringnum': [2, 3, 4], 'type': ['Si', 'O']} 
    
    id =  8  ptclObj =   O 8.000000 2.300000 -20.100000 1.000000 8.000000    molnum  2 ringnum 3
    id =  1  ptclObj =   Si 0.200000 1.300000 33.000000 2.000000 1.230000    molnum  1 ringnum 4
    id =  4  ptclObj =   Si 0.000000 2.300000 -20.100000 1.000000 2.340000    molnum  2 ringnum 4
    id =  6  ptclObj =   Si 8.000000 2.300000 -20.100000 1.000000 8.000000    molnum  2 ringnum 3
    id =  7  ptclObj =   O 8.000000 2.300000 -20.100000 1.000000 8.000000    molnum  2 ringnum 3

