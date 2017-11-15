.. _forcefields_example:
  
forcefields_example
========================
 

This notebook imports the fundamental objects of the streamm.forcefields
module and goes through the functionality of each

.. code:: python

    from pprint import pprint

.. code:: python

    import logging
    logging.basicConfig(filename='forcefield_example.log',level=logging.DEBUG)

.. code:: python

    from pathlib2 import Path
    import os

.. code:: python

    import streamm.forcefields.particletype as particletype

Let’s start with an ethane molecule from the
buildingblocks_example.ipynb example

We have an sp3 carbon bonded to hydrogens

Let’s create the force field parameters object for a ‘CT’ carbon and an
‘HC’ hydrogen

.. code:: python

    CT = particletype.Particletype('CT')
    HC = particletype.Particletype('HC')

Set some parameters from J. Am. Chem. Soc., 1996, 118 (45), pp
11225–11236

http://pubs.acs.org/doi/suppl/10.1021/ja9621760/suppl_file/ja11225.pdf

In general you should pick a force field that has been shown to work
well for your system and set up the parameters

Check that we have our units set right

.. code:: python

    print CT.unit_conf['energy'],CT.unit_conf['length']

Our potential is in kCal/mol (``kCalmol``) so let’s get the unit
dictionary and create our own default

.. code:: python

    import copy

.. code:: python

    oplsaa_unit_conf = copy.deepcopy(CT.unit_conf)

.. code:: python

    oplsaa_unit_conf['energy'] = 'kCalmol'

.. code:: python

    pprint(oplsaa_unit_conf)

.. code:: python

    CT.update_units(oplsaa_unit_conf)

.. code:: python

    HC.update_units(oplsaa_unit_conf)

.. code:: python

    CT.epsilon = 0.066 # kcal/mol
    CT.sigma = 3.5 # Angstroms 

.. code:: python

    HC.epsilon = 0.03 # kcal/mol
    HC.sigma = 2.5 # Angstroms 

Set mass using periodic table

.. code:: python

    import pymatgen_core.core.periodic_table as periodic_table

.. code:: python

    CT.mass =  periodic_table.Element['C'].atomic_mass.real
    HC.mass =  periodic_table.Element['H'].atomic_mass.real

Set the bond stretching parameters

.. code:: python

    import streamm.forcefields.bondtype as bondtype

.. code:: python

    C_H = bondtype.Bondtype('CT','HC',unit_conf=oplsaa_unit_conf) 

.. code:: python

    C_H.setharmonic(1.080,367.0)

.. code:: python

    print C_H

.. code:: python

    C_C = bondtype.Bondtype('CT','CT',unit_conf=oplsaa_unit_conf)
    C_C.setharmonic(1.53,268.0)

.. code:: python

    import streamm.forcefields.angletype as angletype

.. code:: python

    H_C_H = angletype.Angletype('HC','CT','HC',unit_conf=oplsaa_unit_conf)

.. code:: python

    H_C_H.setharmonic(110.7,37.50)

.. code:: python

    print H_C_H

.. code:: python

    H_C_C = angletype.Angletype('HC','CT','CT',unit_conf=oplsaa_unit_conf)
    H_C_C.setharmonic(110.7,37.50)

Now we need a dihedral potential for the HC-CT-CT-HC dihedral

.. code:: python

    import streamm.forcefields.dihtype as dihtype

.. code:: python

    H_C_C_H = dihtype.Dihtype('HC','CT','CT','HC',unit_conf=oplsaa_unit_conf)

.. code:: python

    H_C_C_H.type ='opls'

.. code:: python

    H_C_C_H.setopls(0.0,0.0,0.3,0.0)

Let’s create a parameter container to keep track of our parameters

.. code:: python

    import streamm.forcefields.parameters as parameters 

.. code:: python

    paramC = parameters.Parameters('oplsaa',unit_conf=oplsaa_unit_conf)

Add parameters to the container

.. code:: python

    paramC.add_particletype(CT)

.. code:: python

    paramC.add_particletype(HC)

.. code:: python

    paramC.add_bondtype(C_H)
    paramC.add_bondtype(C_C)

.. code:: python

    paramC.add_angletype(H_C_H)
    paramC.add_angletype(H_C_C)

.. code:: python

    paramC.add_dihtype(H_C_C_H)

.. code:: python

    print paramC

.. code:: python

    for ptkey,pt in paramC.particletypes.iteritems():
        print ptkey,pt,pt.unit_conf['energy'],pt.unit_conf['length']
        

.. code:: python

    for btkey,bt in paramC.bondtypes.iteritems():
        print btkey,bt,bt.unit_conf['harm_bond_coeff'],pt.unit_conf['length']

.. code:: python

    for atkey,at in paramC.angletypes.iteritems():
        print atkey,at,at.unit_conf['energy'],at.unit_conf['length']

.. code:: python

    print paramC.tag

.. code:: python

    paramC.unit_conf

We can dump a pickle file

.. code:: python

    paramC.dump_pickle()

Or we can export a json object

.. code:: python

    paramC_json = paramC.export_json()

Read in ethane .json file from the structures example

.. code:: python

    import streamm.structures.buildingblock as bb

.. code:: python

    need_files = ['ethane_struc.json']
    for f in need_files:
        path = Path(f)
        if not path.is_file():
            print("Need to run buildingblocks_example.ipynb")
            os.system("jupyter nbconvert --to python  buildingblocks_example.ipynb")
            os.system("python buildingblocks_example.py")

.. code:: python

    mol_i = bb.Buildingblock('ethane')

.. code:: python

    mol_i.import_json()

.. code:: python

    print(mol_i.print_properties())

Let’s set the ``paramkey`` for each particle based on the symbol.

.. code:: python

    for pk,p in mol_i.particles.iteritems():
        print  p.symbol 
        if( p.symbol == 'C' ):
            p.paramkey = 'CA'
        elif( p.symbol == 'H' ):
            p.paramkey = 'HA' 
        print p.paramkey ,mol_i.bonded_nblist.calc_nnab(pk)


This is a bit redundant, but we can think of a more complex molecule
where we could use the number of neighbors to write a more complex
routine

.. code:: python

    print mol_i.n_particles

Now we can set the particles, bonds, bond angles and dihedrals of the
molecule to have parameters

First lets set the particle types

.. code:: python

    for pk,p in mol_i.particles.iteritems():
        if( p.paramkey == 'CA' ):
            p.param = CT
            p.param_index = 0
        elif( p.paramkey == 'HA' ):
            p.param = HC
            p.param_index = 1


Now we can set the bond types

.. code:: python

    for bk,b in mol_i.bonds.iteritems():
        b.param = C_H
        b.param_index = 0 

.. code:: python

    for ak,a in mol_i.angles.iteritems():
        a.param = H_C_H
        b.param_index = 0 

.. code:: python

    for dk,d in mol_i.dihedrals.iteritems():
        d.param = H_C_C_H
        d.param_index = 0 

.. code:: python

    print "Particles "
    for pk,p in mol_i.particles.iteritems():
        print p,p.param, p.param_index 
    print "\n Bonds "
    for bk,b in mol_i.bonds.iteritems():    
        print b,b.param, b.param_index 
    print "\n Bond angles "
    for ak,a in mol_i.angles.iteritems():
        print a,a.param, a.param_index 
    print "\n Dihedrals "
    for ak,a in mol_i.dihedrals.iteritems():
        print a,a.param, a.param_index     

Now our molecule has forcefield parameters for all the interactions

Now let’s say we want to use a software package like GROMACS that used
kJ/mol instead of kCal/mol

.. code:: python

    gromacs_unit_conf = copy.deepcopy(oplsaa_unit_conf)

.. code:: python

    gromacs_unit_conf['energy'] = 'kJmol'
    gromacs_unit_conf['length'] = 'nm'
    
    gromacs_unit_conf['harm_bond_coeff'] = 'kJmolsqnm' #*

-  The harmonic bond coefficient ``harm_bond_coeff`` has to be changed
   as well since it has special units of energy/length^2

.. code:: python

    pprint(gromacs_unit_conf)

.. code:: python

    mol_i.update_units(gromacs_unit_conf)

.. code:: python

    print "Particles "
    for pk,p in mol_i.particles.iteritems():
        print p,p.param, p.param_index 
    print "\n Bonds "
    for bk,b in mol_i.bonds.iteritems():    
        print b,b.param, b.param_index 
    print "\n Bond angles "
    for ak,a in mol_i.angles.iteritems():
        print a,a.param, a.param_index 
    print "\n Dihedrals "
    for ak,a in mol_i.dihedrals.iteritems():
        print a,a.param, a.param_index      

.. code:: python

    mol_i.update_units(oplsaa_unit_conf)

.. code:: python

    print "Particles "
    for pk,p in mol_i.particles.iteritems():
        print p,p.param, p.param_index 
    print "\n Bonds "
    for bk,b in mol_i.bonds.iteritems():    
        print b,b.param, b.param_index 
    print "\n Bond angles "
    for ak,a in mol_i.angles.iteritems():
        print a,a.param, a.param_index 
    print "\n Dihedrals "
    for ak,a in mol_i.dihedrals.iteritems():
        print a,a.param, a.param_index      

Sweet as, bro!
