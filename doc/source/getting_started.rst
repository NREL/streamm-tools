.. _getting_started:

Getting started
***************

The streamm package includes modules to manipulate atomic structures, keep track of force field parameters and write input files for Quantum chemistry and molecular dynamics codes.

First create a streamm Builingblock of methane from scratch. ::

    import streamm
    methane = streamm.Buildingblock('methane')
    C = streamm.Particle(symbol='C')
    H = streamm.Particle(symbol='H')
    methane.add_partpos(C,[0.0,0.0,0.0])
    methane.add_partpos(H,[0.69,0.69,0.69])
    methane.add_partpos(H,[-0.69,-0.69,0.69])
    methane.add_partpos(H,[-0.69,0.69,-0.69])
    methane.add_partpos(H,[0.69,-0.69,-0.69])

You could also us a molecular viewer such as `Avogadro <https://avogadro.cc/>`_ and export an .xyz file,
then read the .xyz file into a streamm Builingblock object.
Next build a neighbor list based on the bonded_radius of each particle,
and find all the bonds, bond angles and dihedrals of the Buildingblock. ::
    
    methane.bonded_nblist = methane.guess_nblist(0,radii_buffer=1.25)
    methane.bonded_bonds()
    methane.bonded_angles()
    methane.bonded_dih()
    
Now label some hydrogens as particles that can be subsituded during a reaction. ::
    
    methane.particles[1].rsite = 'RH'
    methane.particles[2].rsite = 'RH'
    methane.find_rsites()

We labeled these sites as 'RH', but it does not really matter, as long as you pass these labels to the attach function. ::

    import streamm.structures.buildingblock as bb
    ethane = bb.attach(methane,methane,'RH',0,'RH',1,tag='ethane')

If we want to run some MD using force fields we need to set up a Parameters container. ::

    oplsaa = streamm.Parameters('oplsaa')

Let's set the energy and legnth units we will input from the literature.::

    oplsaa.update_units({'energy':'kCalmol','length':'ang'})
    
Add some particle types to our Parameters container and pass in the units we are using. ::
    
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

Add some bond, bond angle and dihedral types. ::
    
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

Now we need to set the paramkey of each particle in are Buildingblock to have a key matching a Particletype key. ::
    
    for pk,p in ethane.particles.iteritems():
        if( p.symbol == 'C' ):
            p.paramkey = 'CT'
        elif( p.symbol == 'H' ):
            p.paramkey = 'HC' 

If want to run a LAMMPS simulation we can create a calculation object. ::

    md_calc = streamm.LAMMPS('ethane_md')
    
Set our Buildingblock and Parameter objects to have the correct units for a LAMMPS simulation
and add the Calculation object. ::
    
    ethane.update_units(md_calc.unit_conf)
    oplsaa.update_units(md_calc.unit_conf)
    md_calc.strucC = ethane
    md_calc.paramC = oplsaa

Then we can use the set_ffparam() function to match all the force field parameters to the Buildingblock
based on their paramkeys. ::

    md_calc.set_ffparam()
        
Finally we can output a LAMMPS .data input file for our calculation. ::

    md_calc.write_data()
    



