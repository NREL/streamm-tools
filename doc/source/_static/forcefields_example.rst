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

In general, you should pick a force field that has been shown to work
well for your system and set up the parameters

Check that we have our units set right

.. code:: python

    print CT.unit_conf['energy'],CT.unit_conf['length']


.. parsed-literal::

    Ha ang


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


.. parsed-literal::

    {u'amount': u'atom',
     u'angle': u'degree',
     u'capacitance': u'F',
     u'charge': u'e',
     u'conductance': u'S',
     u'current': u'A',
     u'density': u'amu_nm^3',
     u'electric_dipole_moment': u'D',
     u'emf': u'V',
     u'energy': 'kCalmol',
     u'force': u'GN',
     u'frequency': u'Hz',
     u'harm_bond_coeff': u'kCalmolsqang',
     u'intensity': u'cd',
     u'length': u'ang',
     u'magnetic_flux': u'Wb',
     u'mass': u'amu',
     u'memory': u'Kb',
     u'power': u'GW',
     u'pressure': u'KPa',
     u'resistance': u'ohm',
     u'temperature': u'K',
     u'time': u'ns',
     u'volume': u'nm^3'}


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


.. parsed-literal::

     bond  CT - HC type harmonic 
      harmonic r_0 = 1.080000 K = 367.000000 lammps index 0  gromacs index 0  


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


.. parsed-literal::

     angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0  


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


.. parsed-literal::

    
        Parameters 
          LJ parameters 2 
          Bond parameters 2 
          Angle parameters 2 
          Dihedral parameters 1 
          Improper Dihedral parameters 0 
    


.. code:: python

    for ptkey,pt in paramC.particletypes.iteritems():
        print ptkey,pt,pt.unit_conf['energy'],pt.unit_conf['length']
        


.. parsed-literal::

    0  CT epsilon:0.066 sigma:3.5 kCalmol ang
    1  HC epsilon:0.03 sigma:2.5 kCalmol ang


.. code:: python

    for btkey,bt in paramC.bondtypes.iteritems():
        print btkey,bt,bt.unit_conf['harm_bond_coeff'],pt.unit_conf['length']


.. parsed-literal::

    0  bond  CT - HC type harmonic 
      harmonic r_0 = 1.080000 K = 367.000000 lammps index 0  gromacs index 0   kCalmolsqang ang
    1  bond  CT - CT type harmonic 
      harmonic r_0 = 1.530000 K = 268.000000 lammps index 0  gromacs index 0   kCalmolsqang ang


.. code:: python

    for atkey,at in paramC.angletypes.iteritems():
        print atkey,at,at.unit_conf['energy'],at.unit_conf['length']


.. parsed-literal::

    0  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   kCalmol ang
    1  angle  HC - CT - CT type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   kCalmol ang


.. code:: python

    print paramC.tag


.. parsed-literal::

    oplsaa


.. code:: python

    paramC.unit_conf




.. parsed-literal::

    {u'amount': u'atom',
     u'angle': u'degree',
     u'capacitance': u'F',
     u'charge': u'e',
     u'conductance': u'S',
     u'current': u'A',
     u'density': u'amu_nm^3',
     u'electric_dipole_moment': u'D',
     u'emf': u'V',
     u'energy': 'kCalmol',
     u'force': u'GN',
     u'frequency': u'Hz',
     u'harm_bond_coeff': u'kCalmolsqang',
     u'intensity': u'cd',
     u'length': u'ang',
     u'magnetic_flux': u'Wb',
     u'mass': u'amu',
     u'memory': u'Kb',
     u'power': u'GW',
     u'pressure': u'KPa',
     u'resistance': u'ohm',
     u'temperature': u'K',
     u'time': u'ns',
     u'volume': u'nm^3'}



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


.. parsed-literal::

     n_particles:8 
     n_bonds:7
     n_angles:12
     n_dihedrals:9
     n_impropers:0


Let’s set the ``paramkey`` for each particle based on the symbol.

.. code:: python

    for pk,p in mol_i.particles.iteritems():
        print  p.symbol 
        if( p.symbol == 'C' ):
            p.paramkey = 'CA'
        elif( p.symbol == 'H' ):
            p.paramkey = 'HA' 
        print p.paramkey ,mol_i.bonded_nblist.calc_nnab(pk)



.. parsed-literal::

    C
    CA 4
    H
    HA 1
    H
    HA 1
    H
    HA 1
    C
    CA 4
    H
    HA 1
    H
    HA 1
    H
    HA 1


This is a bit redundant, but we can think of a more complex molecule
where we could use the number of neighbors to write a more complex
routine

.. code:: python

    print mol_i.n_particles


.. parsed-literal::

    8


Now we can set the particles, bonds, bond angles and dihedrals of the
molecule to have parameters

First let's set the particle types

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


.. parsed-literal::

    Particles 
    atom C (C)  CT epsilon:0.066 sigma:3.5 0
    atom H (H)  HC epsilon:0.03 sigma:2.5 1
    atom H (H)  HC epsilon:0.03 sigma:2.5 1
    atom H (H)  HC epsilon:0.03 sigma:2.5 1
    atom C (C)  CT epsilon:0.066 sigma:3.5 0
    atom H (H)  HC epsilon:0.03 sigma:2.5 1
    atom H (H)  HC epsilon:0.03 sigma:2.5 1
    atom H (H)  HC epsilon:0.03 sigma:2.5 1
    
     Bonds 
     0 - 1  bond  CT - HC type harmonic 
      harmonic r_0 = 1.080000 K = 367.000000 lammps index 0  gromacs index 0   0
     0 - 2  bond  CT - HC type harmonic 
      harmonic r_0 = 1.080000 K = 367.000000 lammps index 0  gromacs index 0   0
     0 - 3  bond  CT - HC type harmonic 
      harmonic r_0 = 1.080000 K = 367.000000 lammps index 0  gromacs index 0   0
     0 - 4  bond  CT - HC type harmonic 
      harmonic r_0 = 1.080000 K = 367.000000 lammps index 0  gromacs index 0   0
     4 - 5  bond  CT - HC type harmonic 
      harmonic r_0 = 1.080000 K = 367.000000 lammps index 0  gromacs index 0   0
     4 - 6  bond  CT - HC type harmonic 
      harmonic r_0 = 1.080000 K = 367.000000 lammps index 0  gromacs index 0   0
     4 - 7  bond  CT - HC type harmonic 
      harmonic r_0 = 1.080000 K = 367.000000 lammps index 0  gromacs index 0   0
    
     Bond angles 
     2 - 0 - 1  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   0
     3 - 0 - 1  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   0
     4 - 0 - 1  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   0
     3 - 0 - 2  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   0
     4 - 0 - 2  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   0
     4 - 0 - 3  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   0
     5 - 4 - 0  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   0
     6 - 4 - 0  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   0
     7 - 4 - 0  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   0
     6 - 4 - 5  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   0
     7 - 4 - 5  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   0
     7 - 4 - 6  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   0
    
     Dihedrals 
     1 - 0 - 4 - 5  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 0.300000 k4 = 0.000000 lammps index 0  gromcas index 0  0
     1 - 0 - 4 - 6  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 0.300000 k4 = 0.000000 lammps index 0  gromcas index 0  0
     1 - 0 - 4 - 7  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 0.300000 k4 = 0.000000 lammps index 0  gromcas index 0  0
     2 - 0 - 4 - 5  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 0.300000 k4 = 0.000000 lammps index 0  gromcas index 0  0
     2 - 0 - 4 - 6  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 0.300000 k4 = 0.000000 lammps index 0  gromcas index 0  0
     2 - 0 - 4 - 7  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 0.300000 k4 = 0.000000 lammps index 0  gromcas index 0  0
     3 - 0 - 4 - 5  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 0.300000 k4 = 0.000000 lammps index 0  gromcas index 0  0
     3 - 0 - 4 - 6  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 0.300000 k4 = 0.000000 lammps index 0  gromcas index 0  0
     3 - 0 - 4 - 7  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 0.300000 k4 = 0.000000 lammps index 0  gromcas index 0  0


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


.. parsed-literal::

    {u'amount': u'atom',
     u'angle': u'degree',
     u'capacitance': u'F',
     u'charge': u'e',
     u'conductance': u'S',
     u'current': u'A',
     u'density': u'amu_nm^3',
     u'electric_dipole_moment': u'D',
     u'emf': u'V',
     u'energy': 'kJmol',
     u'force': u'GN',
     u'frequency': u'Hz',
     u'harm_bond_coeff': 'kJmolsqnm',
     u'intensity': u'cd',
     u'length': 'nm',
     u'magnetic_flux': u'Wb',
     u'mass': u'amu',
     u'memory': u'Kb',
     u'power': u'GW',
     u'pressure': u'KPa',
     u'resistance': u'ohm',
     u'temperature': u'K',
     u'time': u'ns',
     u'volume': u'nm^3'}


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


.. parsed-literal::

    Particles 
    atom C (C)  CT epsilon:0.276144 sigma:0.35 0
    atom H (H)  HC epsilon:0.12552 sigma:0.25 1
    atom H (H)  HC epsilon:0.12552 sigma:0.25 1
    atom H (H)  HC epsilon:0.12552 sigma:0.25 1
    atom C (C)  CT epsilon:0.276144 sigma:0.35 0
    atom H (H)  HC epsilon:0.12552 sigma:0.25 1
    atom H (H)  HC epsilon:0.12552 sigma:0.25 1
    atom H (H)  HC epsilon:0.12552 sigma:0.25 1
    
     Bonds 
     0 - 1  bond  CT - HC type harmonic 
      harmonic r_0 = 0.108000 K = 153552.800000 lammps index 0  gromacs index 0   0
     0 - 2  bond  CT - HC type harmonic 
      harmonic r_0 = 0.108000 K = 153552.800000 lammps index 0  gromacs index 0   0
     0 - 3  bond  CT - HC type harmonic 
      harmonic r_0 = 0.108000 K = 153552.800000 lammps index 0  gromacs index 0   0
     0 - 4  bond  CT - HC type harmonic 
      harmonic r_0 = 0.108000 K = 153552.800000 lammps index 0  gromacs index 0   0
     4 - 5  bond  CT - HC type harmonic 
      harmonic r_0 = 0.108000 K = 153552.800000 lammps index 0  gromacs index 0   0
     4 - 6  bond  CT - HC type harmonic 
      harmonic r_0 = 0.108000 K = 153552.800000 lammps index 0  gromacs index 0   0
     4 - 7  bond  CT - HC type harmonic 
      harmonic r_0 = 0.108000 K = 153552.800000 lammps index 0  gromacs index 0   0
    
     Bond angles 
     2 - 0 - 1  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 156.900000 lammps index 0  gromacs index 0   0
     3 - 0 - 1  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 156.900000 lammps index 0  gromacs index 0   0
     4 - 0 - 1  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 156.900000 lammps index 0  gromacs index 0   0
     3 - 0 - 2  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 156.900000 lammps index 0  gromacs index 0   0
     4 - 0 - 2  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 156.900000 lammps index 0  gromacs index 0   0
     4 - 0 - 3  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 156.900000 lammps index 0  gromacs index 0   0
     5 - 4 - 0  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 156.900000 lammps index 0  gromacs index 0   0
     6 - 4 - 0  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 156.900000 lammps index 0  gromacs index 0   0
     7 - 4 - 0  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 156.900000 lammps index 0  gromacs index 0   0
     6 - 4 - 5  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 156.900000 lammps index 0  gromacs index 0   0
     7 - 4 - 5  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 156.900000 lammps index 0  gromacs index 0   0
     7 - 4 - 6  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 156.900000 lammps index 0  gromacs index 0   0
    
     Dihedrals 
     1 - 0 - 4 - 5  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 1.255200 k4 = 0.000000 lammps index 0  gromcas index 0  0
     1 - 0 - 4 - 6  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 1.255200 k4 = 0.000000 lammps index 0  gromcas index 0  0
     1 - 0 - 4 - 7  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 1.255200 k4 = 0.000000 lammps index 0  gromcas index 0  0
     2 - 0 - 4 - 5  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 1.255200 k4 = 0.000000 lammps index 0  gromcas index 0  0
     2 - 0 - 4 - 6  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 1.255200 k4 = 0.000000 lammps index 0  gromcas index 0  0
     2 - 0 - 4 - 7  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 1.255200 k4 = 0.000000 lammps index 0  gromcas index 0  0
     3 - 0 - 4 - 5  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 1.255200 k4 = 0.000000 lammps index 0  gromcas index 0  0
     3 - 0 - 4 - 6  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 1.255200 k4 = 0.000000 lammps index 0  gromcas index 0  0
     3 - 0 - 4 - 7  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 1.255200 k4 = 0.000000 lammps index 0  gromcas index 0  0


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


.. parsed-literal::

    Particles 
    atom C (C)  CT epsilon:0.066 sigma:3.5 0
    atom H (H)  HC epsilon:0.03 sigma:2.5 1
    atom H (H)  HC epsilon:0.03 sigma:2.5 1
    atom H (H)  HC epsilon:0.03 sigma:2.5 1
    atom C (C)  CT epsilon:0.066 sigma:3.5 0
    atom H (H)  HC epsilon:0.03 sigma:2.5 1
    atom H (H)  HC epsilon:0.03 sigma:2.5 1
    atom H (H)  HC epsilon:0.03 sigma:2.5 1
    
     Bonds 
     0 - 1  bond  CT - HC type harmonic 
      harmonic r_0 = 1.080000 K = 367.000000 lammps index 0  gromacs index 0   0
     0 - 2  bond  CT - HC type harmonic 
      harmonic r_0 = 1.080000 K = 367.000000 lammps index 0  gromacs index 0   0
     0 - 3  bond  CT - HC type harmonic 
      harmonic r_0 = 1.080000 K = 367.000000 lammps index 0  gromacs index 0   0
     0 - 4  bond  CT - HC type harmonic 
      harmonic r_0 = 1.080000 K = 367.000000 lammps index 0  gromacs index 0   0
     4 - 5  bond  CT - HC type harmonic 
      harmonic r_0 = 1.080000 K = 367.000000 lammps index 0  gromacs index 0   0
     4 - 6  bond  CT - HC type harmonic 
      harmonic r_0 = 1.080000 K = 367.000000 lammps index 0  gromacs index 0   0
     4 - 7  bond  CT - HC type harmonic 
      harmonic r_0 = 1.080000 K = 367.000000 lammps index 0  gromacs index 0   0
    
     Bond angles 
     2 - 0 - 1  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   0
     3 - 0 - 1  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   0
     4 - 0 - 1  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   0
     3 - 0 - 2  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   0
     4 - 0 - 2  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   0
     4 - 0 - 3  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   0
     5 - 4 - 0  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   0
     6 - 4 - 0  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   0
     7 - 4 - 0  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   0
     6 - 4 - 5  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   0
     7 - 4 - 5  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   0
     7 - 4 - 6  angle  HC - CT - HC type harmonic 
      harmonic theta_0 = 110.700000 K = 37.500000 lammps index 0  gromacs index 0   0
    
     Dihedrals 
     1 - 0 - 4 - 5  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 0.300000 k4 = 0.000000 lammps index 0  gromcas index 0  0
     1 - 0 - 4 - 6  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 0.300000 k4 = 0.000000 lammps index 0  gromcas index 0  0
     1 - 0 - 4 - 7  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 0.300000 k4 = 0.000000 lammps index 0  gromcas index 0  0
     2 - 0 - 4 - 5  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 0.300000 k4 = 0.000000 lammps index 0  gromcas index 0  0
     2 - 0 - 4 - 6  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 0.300000 k4 = 0.000000 lammps index 0  gromcas index 0  0
     2 - 0 - 4 - 7  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 0.300000 k4 = 0.000000 lammps index 0  gromcas index 0  0
     3 - 0 - 4 - 5  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 0.300000 k4 = 0.000000 lammps index 0  gromcas index 0  0
     3 - 0 - 4 - 6  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 0.300000 k4 = 0.000000 lammps index 0  gromcas index 0  0
     3 - 0 - 4 - 7  dihedral  HC - CT - CT - HC type opls 
      k1 = 0.000000 k2 = 0.000000 k3 = 0.300000 k4 = 0.000000 lammps index 0  gromcas index 0  0


Sweet as, bro!
