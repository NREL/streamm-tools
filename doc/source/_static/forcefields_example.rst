.. _forcefields_example:
  
forcefields_example
===============
 

This notebook imports the fundamental objects of the streamm.forcefields
module and goes through the functionality of each

.. code:: ipython2

    from pprint import pprint

.. code:: ipython2

    import logging
    logging.basicConfig(filename='forcefield_example.log',level=logging.DEBUG)

.. code:: ipython2

    from pathlib2 import Path
    import os

.. code:: ipython2

    import streamm.forcefields.particletype as particletype

Let’s start with an ethane molecule from the
buildingblocks_example.ipynb example

We have an sp3 carbon bonded to hydrogens

Let’s create the force field parameters object for a ‘CT’ carbon and an
‘HC’ hydrogen

.. code:: ipython2

    CT = particletype.Particletype('CT')
    HC = particletype.Particletype('HC')

Set some parameters from J. Am. Chem. Soc., 1996, 118 (45), pp
11225–11236

http://pubs.acs.org/doi/suppl/10.1021/ja9621760/suppl_file/ja11225.pdf

In general you should pick a force field that has been shown to work
well for your system and set up the parameters

Check that we have our units set right

.. code:: ipython2

    print CT.unit_conf['energy'],CT.unit_conf['length']

Our potential is in kCal/mol (``kCalmol``) so let’s get the unit
dictionary and create our own default

.. code:: ipython2

    import copy

.. code:: ipython2

    oplsaa_unit_conf = copy.deepcopy(CT.unit_conf)

.. code:: ipython2

    oplsaa_unit_conf['energy'] = 'kCalmol'

.. code:: ipython2

    pprint(oplsaa_unit_conf)

.. code:: ipython2

    CT.update_units(oplsaa_unit_conf)

.. code:: ipython2

    HC.update_units(oplsaa_unit_conf)

.. code:: ipython2

    CT.epsilon = 0.066 # kcal/mol
    CT.sigma = 3.5 # Angstroms 

.. code:: ipython2

    HC.epsilon = 0.03 # kcal/mol
    HC.sigma = 2.5 # Angstroms 

Set mass using periodic table

.. code:: ipython2

    import pymatgen_core.core.periodic_table as periodic_table

.. code:: ipython2

    CT.mass =  periodic_table.Element['C'].atomic_mass.real
    HC.mass =  periodic_table.Element['H'].atomic_mass.real

Set the bond stretching parameters

.. code:: ipython2

    import streamm.forcefields.bondtype as bondtype

.. code:: ipython2

    C_H = bondtype.Bondtype('CT','HC',unit_conf=oplsaa_unit_conf) 

.. code:: ipython2

    C_H.setharmonic(1.080,367.0)

.. code:: ipython2

    print C_H

.. code:: ipython2

    C_C = bondtype.Bondtype('CT','CT',unit_conf=oplsaa_unit_conf)
    C_C.setharmonic(1.53,268.0)

.. code:: ipython2

    import streamm.forcefields.angletype as angletype

.. code:: ipython2

    H_C_H = angletype.Angletype('HC','CT','HC',unit_conf=oplsaa_unit_conf)

.. code:: ipython2

    H_C_H.setharmonic(110.7,37.50)

.. code:: ipython2

    print H_C_H

.. code:: ipython2

    H_C_C = angletype.Angletype('HC','CT','CT',unit_conf=oplsaa_unit_conf)
    H_C_C.setharmonic(110.7,37.50)

Now we need a dihedral potential for the HC-CT-CT-HC dihedral

.. code:: ipython2

    import streamm.forcefields.dihtype as dihtype

.. code:: ipython2

    H_C_C_H = dihtype.Dihtype('HC','CT','CT','HC',unit_conf=oplsaa_unit_conf)

.. code:: ipython2

    H_C_C_H.type ='opls'

.. code:: ipython2

    H_C_C_H.setopls(0.0,0.0,0.3,0.0)

Let’s create a parameter container to keep track of our parameters

.. code:: ipython2

    import streamm.forcefields.parameters as parameters 

.. code:: ipython2

    paramC = parameters.Parameters('oplsaa',unit_conf=oplsaa_unit_conf)

Add parameters to the container

.. code:: ipython2

    paramC.add_particletype(CT)

.. code:: ipython2

    paramC.add_particletype(HC)

.. code:: ipython2

    paramC.add_bondtype(C_H)
    paramC.add_bondtype(C_C)

.. code:: ipython2

    paramC.add_angletype(H_C_H)
    paramC.add_angletype(H_C_C)

.. code:: ipython2

    paramC.add_dihtype(H_C_C_H)

.. code:: ipython2

    print paramC

.. code:: ipython2

    for ptkey,pt in paramC.particletypes.iteritems():
        print ptkey,pt,pt.unit_conf['energy'],pt.unit_conf['length']
        

.. code:: ipython2

    for btkey,bt in paramC.bondtypes.iteritems():
        print btkey,bt,bt.unit_conf['harm_bond_coeff'],pt.unit_conf['length']

.. code:: ipython2

    for atkey,at in paramC.angletypes.iteritems():
        print atkey,at,at.unit_conf['energy'],at.unit_conf['length']

.. code:: ipython2

    print paramC.tag

.. code:: ipython2

    paramC.unit_conf

We can dump a pickle file

.. code:: ipython2

    paramC.dump_pickle()

Or we can export a json object

.. code:: ipython2

    paramC_json = paramC.export_json()

Read in ethane .json file from the structures example

.. code:: ipython2

    import streamm.structures.buildingblock as bb

.. code:: ipython2

    need_files = ['ethane_struc.json']
    for f in need_files:
        path = Path(f)
        if not path.is_file():
            print("Need to run buildingblocks_example.ipynb")
            os.system("jupyter nbconvert --to python  buildingblocks_example.ipynb")
            os.system("python buildingblocks_example.py")

.. code:: ipython2

    mol_i = bb.Buildingblock('ethane')

.. code:: ipython2

    mol_i.import_json()

.. code:: ipython2

    print(mol_i.print_properties())

Let’s set the ``paramkey`` for each particle based on the symbol.

.. code:: ipython2

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

.. code:: ipython2

    print mol_i.n_particles

Now we can set the particles, bonds, bond angles and dihedrals of the
molecule to have parameters

First lets set the particle types

.. code:: ipython2

    for pk,p in mol_i.particles.iteritems():
        if( p.paramkey == 'CA' ):
            p.param = CT
            p.param_index = 0
        elif( p.paramkey == 'HA' ):
            p.param = HC
            p.param_index = 1


Now we can set the bond types

.. code:: ipython2

    for bk,b in mol_i.bonds.iteritems():
        b.param = C_H
        b.param_index = 0 

.. code:: ipython2

    for ak,a in mol_i.angles.iteritems():
        a.param = H_C_H
        b.param_index = 0 

.. code:: ipython2

    for dk,d in mol_i.dihedrals.iteritems():
        d.param = H_C_C_H
        d.param_index = 0 

.. code:: ipython2

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

.. code:: ipython2

    gromacs_unit_conf = copy.deepcopy(oplsaa_unit_conf)

.. code:: ipython2

    gromacs_unit_conf['energy'] = 'kJmol'
    gromacs_unit_conf['length'] = 'nm'
    
    gromacs_unit_conf['harm_bond_coeff'] = 'kJmolsqnm' #*

-  The harmonic bond coefficient ``harm_bond_coeff`` has to be changed
   as well since it has special units of energy/length^2

.. code:: ipython2

    pprint(gromacs_unit_conf)

.. code:: ipython2

    mol_i.update_units(gromacs_unit_conf)

.. code:: ipython2

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

.. code:: ipython2

    mol_i.update_units(oplsaa_unit_conf)

.. code:: ipython2

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
