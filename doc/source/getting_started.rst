.. _getting_started:

Getting started
***************

The streamm package includes modules to manipulate atomic structures,
keep track of force field parameters and write input
files for Quantum chemistry codes, such as
`NWChem <http://www.nwchem-sw.org/index.php/Main_Page>`_
and `Gaussian <http://gaussian.com/>`_ molecular dynamics codes
such as `LAMMPS <http://lammps.sandia.gov/>`_.

To get started we will first build ethane out of two methanes using the :class:`streamm.Buildingblock <streamm.structures.buildingblock.Buildingblock>` object.
While this is a simple example, it is meant to illustrate the functionality of the streamm code to construct organic structures form building block units.
In practice, this functionality allows for the combinatorial creation and analysis millions of variations of organic structures.

To create a :class:`streamm.Buildingblock <streamm.structures.buildingblock.Buildingblock>`
object representing methane, we will create carbon and hydrogen particle objects and add them to the methane object with the correct positions.

.. code :: python 

    import streamm
    methane = streamm.Buildingblock('methane')
    C = streamm.Particle(symbol='C')
    H = streamm.Particle(symbol='H')
    methane.add_partpos(C,[0.0,0.0,0.0])
    methane.add_partpos(H,[0.69,0.69,0.69])
    methane.add_partpos(H,[-0.69,-0.69,0.69])
    methane.add_partpos(H,[-0.69,0.69,-0.69])
    methane.add_partpos(H,[0.69,-0.69,-0.69])


You could also use a molecular viewer such as `Avogadro <https://avogadro.cc/>`_ to create an organic structure, see the :ref:`read_xyz` :ref:`how_to` for more information. 

Next, we need to define the connectivity of structure by guessing a
:class:`neighbor list <streamm.structures.nblist.NBlist>` based on the
`bonded_radius` of each :class:`Particle <streamm.structures.particle.Particle>` using the :func:`guess_nblist() <streamm.structures.structure.Structure.guess_nblist>` function. 
    
.. code :: python 
 
    methane.bonded_nblist = methane.guess_nblist(0,radii_buffer=1.25)
    
Then we can label some hydrogens as substitutable sites (rsite), and run the :func:`find_rsites() <streamm.structures.buildingblock.Buildingblock.find_rsites>` function to update the `funcs` list of the
:class:`streamm.Buildingblock <streamm.structures.buildingblock.Buildingblock>` object.

.. code :: python 

    methane.particles[1].rsite = 'RH'
    methane.particles[2].rsite = 'RH'
    methane.find_rsites()

We labeled these sites as `RH`, but it does not really matter, as long as you pass these labels to the :func:`attach() <streamm.structures.buildingblock.attach>` function. 

.. code :: python 

    import streamm.structures.buildingblock as bb
    ethane = bb.attach(methane,methane,'RH',0,'RH',1,tag='ethane')


Then you can write an `.xyz` file to visualize your new :class:`streamm.Buildingblock <streamm.structures.buildingblock.Buildingblock>` using your favorite molecular viewing software.

.. code :: python

    ethane.write_xyz()
    
    
