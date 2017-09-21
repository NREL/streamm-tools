.. _getting_started:

Getting started
***************

.. TODO ::
    
    Add unit_conf
    

The streamm package includes modules to manipulate atomic structures,
keep track of force field parameters and write input
files for Quantum chemistry codes such as
`NWChem <http://www.nwchem-sw.org/index.php/Main_Page>`_
and `Gaussian <http://gaussian.com/>`_ molecular dynamics codes
such as `LAMMPS <http://lammps.sandia.gov/>`_.


First, let us create a
streamm :class:`Buildingblock <streamm.structures.buildingblock.Buildingblock>`
of methane from scratch.

    import streamm
    methane = streamm.Buildingblock('methane')
    C = streamm.Particle(symbol='C')
    H = streamm.Particle(symbol='H')
    methane.add_partpos(C,[0.0,0.0,0.0])
    methane.add_partpos(H,[0.69,0.69,0.69])
    methane.add_partpos(H,[-0.69,-0.69,0.69])
    methane.add_partpos(H,[-0.69,0.69,-0.69])
    methane.add_partpos(H,[0.69,-0.69,-0.69])

You could also use a molecular viewer such as `Avogadro <https://avogadro.cc/>`_
and export a `.xyz <https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/xyz.html>`_ file,
then read the `.xyz <https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/xyz.html>`_
file into a streamm :class:`Buildingblock <streamm.structures.buildingblock.Buildingblock>`  object.
Next, build a neighbor list based on
the `bonded_radius <pymatgen_core.core.units.Particle.bonded_radius>`
of each `Particle <streamm.structures.particle.Particle>`,
and find all the :class:`Bonds <streamm.structures.bond.Bond>`,
:class:`bond angles <streamm.structures.angle.Angle>` and
`dihedrals <streamm.structures.dihedral.Dihedral>` of
the :class:`Buildingblock <streamm.structures.buildingblock.Buildingblock>`.  ::
    
    methane.bonded_nblist = methane.guess_nblist(0,radii_buffer=1.25)
    methane.bonded_bonds()
    methane.bonded_angles()
    methane.bonded_dih()
    
Now we can label some hydrogens as particles as substitutable sites `rsite`. ::

    methane.particles[1].rsite = 'RH'
    methane.particles[2].rsite = 'RH'
    methane.find_rsites()

We labeled these sites as 'RH', but it does not really matter, as long as you pass these labels to the attach function. ::

    import streamm.structures.buildingblock as bb
    ethane = bb.attach(methane,methane,'RH',0,'RH',1,tag='ethane')


MD setup
********

If we want to run some MD using force fields, we need to set up a :class:`Parameters <streamm.forcefields.parameters.Parameters>` container. ::

    oplsaa = streamm.Parameters('oplsaa')

Let's set the energy and legnth units we will input from the literature. ::

    oplsaa.update_units({'energy':'kCalmol','length':'ang'})
    
Add some :class:`Particletype <streamm.forcefields.particletype.Particletype>` objects
to our :class:`Parameters <streamm.forcefields.parameters.Parameters>`
container and pass in the `units_conf` we are using. ::
    
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
:class:`Dihedraltype <streamm.forcefields.dihedraltype.Dihedraltype>` objects. ::
    
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

Now we need to set the `paramkeys` of each particle in
are :class:`Buildingblock <streamm.structures.buildingblock.Buildingblock>`
to have a key matching a :class:`Particletype <streamm.forcefields.particletype.Particletype>` key.    
    for pk,p in ethane.particles.iteritems():
        if( p.symbol == 'C' ):
            p.paramkey = 'CT'
        elif( p.symbol == 'H' ):
            p.paramkey = 'HC' 

If want to run a `LAMMPS <http://lammps.sandia.gov/>` simulation we can create
a :class:`Calculation <streamm.calculations.calculation.Calculation>` object. ::

    md_calc = streamm.LAMMPS('ethane_md')
    
Set our Buildingblock and :class:`Buildingblock <streamm.structures.buildingblock.Buildingblock>`
objects to have the correct units for a `LAMMPS <http://lammps.sandia.gov/>`_
simulation and add the class:`Calculation <streamm.calculations.calculation.Calculation>` object. ::
    
    ethane.update_units(md_calc.unit_conf)
    oplsaa.update_units(md_calc.unit_conf)
    md_calc.strucC = ethane
    md_calc.paramC = oplsaa

Then we can use the :func:`set_ffparam <streamm.calculations.calculation.Calculation.set_ffparam>` function to match all the force field
parameters to the :class:`Buildingblock <streamm.structures.buildingblock.Buildingblock>`  based on their `paramkeys`. ::

    md_calc.set_ffparam()
        
Finally, we can output a LAMMPS `.data` input file for our calculation. ::

    md_calc.write_data()
    



