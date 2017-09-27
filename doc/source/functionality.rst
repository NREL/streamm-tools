.. _functionality:

Functionality
*************

NWChem input file
==================


First, we need to set the location of the templates directory cloned
from http://github.com/NREL/streamm-tools

If you are running this example in the examples directory in the
streamm-tools repo, the TEMPLATE_DIR should look like this

.. code:: python

    EXAMPLE_DIR = os.getcwd()
    print(EXAMPLE_DIR)
    TEMPLATE_DIR =  os.path.join(EXAMPLE_DIR,'..','templates','')
    print(TEMPLATE_DIR)


If not please set the ``TEMPLATE_DIR`` variable to the location of the
templates

We will use the basic nwchem.nw example template

.. code:: python

    temlate_file = 'nwchem.nw'

Create a NWChem calculation object

.. code:: python

    nwchem = streamm.NWChem('ethane_nw_sp')

Set the structure of the calculation to the ethane object.

.. code:: python

    nwchem.strucC = ethane
    

Get the location of the template file

.. code:: python

    template_path =  os.path.join(TEMPLATE_DIR,temlate_file)
    
    print template_path

Read in the template

.. code:: python

    template_line = nwchem.read_lines(template_path)

Set the properties dictionary to contain the information for our
calculation

.. code:: python

    nwchem.properties['basis'] = '6-31g'
    nwchem.properties['method'] = 'UHF'
    nwchem.properties['charge'] = 0
    nwchem.properties['spin_mult'] = 1
    nwchem.properties['task'] = 'SCF '
    nwchem.properties['coord'] = nwchem.strucC.write_coord()

Do a string replace of the dictionary keys to create an input string

.. code:: python

    input_str = nwchem.replace_keys(template_line,nwchem.properties)

.. code:: python

    file_name = '%s.nw'%(nwchem.tag)
    with open(file_name,"w") as F:
        F.write(input_str)


LAMMPS input file
==================

Setting Parameters
------------------

If we want to run some MD using force fields, we need to set up a :class:`Parameters <streamm.forcefields.parameters.Parameters>` container.

.. code :: python 

    oplsaa = streamm.Parameters('oplsaa')

Let's set the energy and length units we will input from the literature.

.. code :: python 

    oplsaa.update_units({'energy':'kCalmol','length':'ang'})
    
Add some :class:`Particletype <streamm.forcefields.particletype.Particletype>` objects
to our :class:`Parameters <streamm.forcefields.parameters.Parameters>`
container and pass in the `units_conf` we are using.

.. code :: python 
    
    CT = streamm.Particletype('CT',unit_conf=oplsaa.unit_conf)
    CT.epsilon = 0.066 # kcal/mol
    CT.sigma = 3.5 # Angstroms 
    CT.mass = 12.0107
    oplsaa.add_particletype(CT)
    HC = streamm.Particletype('HC',unit_conf=oplsaa.unit_conf)
    HC.epsilon = 0.03 # kcal/mol
    HC.sigma = 2.5 # Angstroms 
    HC.mass = 1.00794
    oplsaa.add_particletype(HC)

Add some :class:`Bondtype <streamm.forcefields.bondtype.Bondtype>`,
:class:`Angletype <streamm.forcefields.angletype.Angletype>`, and 
:class:`Dihedraltype <streamm.forcefields.dihedraltype.Dihedraltype>` objects.

.. code :: python 
    
    C_H = streamm.Bondtype('CT','HC',unit_conf=oplsaa.unit_conf)
    C_H.setharmonic(1.080,367.0)
    oplsaa.add_bondtype(C_H)
    
    C_C = streamm.Bondtype('CT','CT',unit_conf=oplsaa.unit_conf)
    C_C.setharmonic(1.080,367.0)
    oplsaa.add_bondtype(C_C)
    
    H_C_H = streamm.Angletype('HC','CT','HC',unit_conf=oplsaa.unit_conf)
    H_C_H.setharmonic(110.7,37.50)
    oplsaa.add_angletype(H_C_H)
    
    H_C_C = streamm.Angletype('HC','CT','CT',unit_conf=oplsaa.unit_conf)
    H_C_C.setharmonic(90.7,60.50)
    oplsaa.add_angletype(H_C_C)

Setting `paramkeys`
-------------------

Now we need to set the `paramkeys` of each particle in
are :class:`Buildingblock <streamm.structures.buildingblock.Buildingblock>`
to have a key matching a :class:`Particletype <streamm.forcefields.particletype.Particletype>` key.

.. code:: python

    for pk,p in ethane.particles.iteritems():
        if( p.symbol == 'C' ):
            p.paramkey = 'CT'
        elif( p.symbol == 'H' ):
            p.paramkey = 'HC' 

Create LAMMPS Calculation
-------------------------------------

If we want to run a `LAMMPS <http://lammps.sandia.gov/>` simulation, we can create
a :class:`Calculation <streamm.calculations.calculation.Calculation>` object. 

.. code:: python

    md_calc = streamm.LAMMPS('ethane_md')
    
Set our Buildingblock and :class:`Buildingblock <streamm.structures.buildingblock.Buildingblock>`
objects to have the correct units for a `LAMMPS <http://lammps.sandia.gov/>`_
simulation and add the class:`Calculation <streamm.calculations.calculation.Calculation>` object.

.. code :: python 
    
    ethane.update_units(md_calc.unit_conf)
    oplsaa.update_units(md_calc.unit_conf)
    md_calc.strucC = ethane
    md_calc.paramC = oplsaa

Find Molecular Connections
------------------------------------

Next, we need to find all the :class:`Bonds <streamm.structures.bond.Bond>`,
:class:`bond angles <streamm.structures.angle.Angle>` and
`dihedrals <streamm.structures.dihedral.Dihedral>` of
the :class:`Buildingblock <streamm.structures.buildingblock.Buildingblock>`, using the bonded :class:`neighbor list <streamm.structures.nblist.NBlist>`.

.. code :: python 
 
     md_calc.strucC.bonded_bonds()
     md_calc.strucC.bonded_angles()
     md_calc.strucC.bonded_dih()

Then we can use the :func:`set_ffparam <streamm.calculations.calculation.Calculation.set_ffparam>` function to match all the force field
parameters to the :class:`Buildingblock <streamm.structures.buildingblock.Buildingblock>`  based on their `paramkeys`.

.. code :: python 

    md_calc.set_ffparam()
        
Finally, we can output a LAMMPS `.data` input file for our calculation.

.. code :: python 

    md_calc.write_data()
    



