
.. code:: python

    import streamm
    
    methane = streamm.Buildingblock('methane')
    C = streamm.Particle(symbol='C')
    H = streamm.Particle(symbol='H')

.. code:: python

    methane.add_partpos(C,[0.0,0.0,0.0])
    methane.add_partpos(H,[0.69,0.69,0.69])
    methane.add_partpos(H,[-0.69,-0.69,0.69])
    methane.add_partpos(H,[-0.69,0.69,-0.69])
    methane.add_partpos(H,[0.69,-0.69,-0.69])

.. code:: python

    methane.bonded_nblist = methane.guess_nblist(0,radii_buffer=1.25)
    
    for p_index,particle_i in methane.particles.iteritems():
        print p_index,particle_i,methane.positions[p_index],methane.bonded_nblist.calc_nnab(p_index),methane.unit_conf['length']
    
    methane.bonded_bonds()
    methane.bonded_angles()
    methane.bonded_dih()
    
    print methane.print_properties()


.. parsed-literal::

    0 atom[0] C (C) [ 0.  0.  0.] 4 ang
    1 atom[1] H (H) [ 0.69  0.69  0.69] 1 ang
    2 atom[2] H (H) [-0.69 -0.69  0.69] 1 ang
    3 atom[3] H (H) [-0.69  0.69 -0.69] 1 ang
    4 atom[4] H (H) [ 0.69 -0.69 -0.69] 1 ang
     n_particles:5 
     n_bonds:4
     n_angles:6
     n_dihedrals:0
     n_impropers:0


.. code:: python

    methane.particles[1].rsite = 'RH'
    methane.particles[2].rsite = 'RH'
    methane.find_rsites()
    print methane.show_rsites()


.. parsed-literal::

    rsite:RH[ paticle:atom[1] H (H) index:1 n_bonds:1] 
    rsite:RH[ paticle:atom[2] H (H) index:2 n_bonds:1] 
    


.. code:: python

    import streamm.structures.buildingblock as bb

.. code:: python

    ethane = bb.attach(methane,methane,'RH',0,'RH',1,tag='ethane')

.. code:: python

    oplsaa = streamm.Parameters('oplsaa')

.. code:: python

    oplsaa.update_units({'energy':'kCalmol','length':'ang'})

.. code:: python

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

.. code:: python

    C_H = streamm.Bondtype('CT','HC',unit_conf=oplsaa.unit_conf)
    C_H.setharmonic(1.08,367.0)
    oplsaa.add_bondtype(C_H)
    
    C_C = streamm.Bondtype('CT','CT',unit_conf=oplsaa.unit_conf)
    C_C.setharmonic(1.53,268.0)
    oplsaa.add_bondtype(C_C)
    
    H_C_H = streamm.Angletype('HC','CT','HC',unit_conf=oplsaa.unit_conf)
    H_C_H.setharmonic(110.7,37.50)
    oplsaa.add_angletype(H_C_H)
    
    H_C_C = streamm.Angletype('HC','CT','CT',unit_conf=oplsaa.unit_conf)
    H_C_C.setharmonic(90.7,60.50)
    oplsaa.add_angletype(H_C_C)

.. code:: python

    for pk,p in ethane.particles.iteritems():
        if( p.symbol == 'C' ):
            p.paramkey = 'CT'
        elif( p.symbol == 'H' ):
            p.paramkey = 'HC' 
        print p.paramkey ,ethane.bonded_nblist.calc_nnab(pk)
    



.. parsed-literal::

    CT 4
    HC 1
    HC 1
    HC 1
    CT 4
    HC 1
    HC 1
    HC 1


.. code:: python

    md_calc = streamm.LAMMPS('ethane_md')

.. code:: python

    ethane.update_units(md_calc.unit_conf)
    
    oplsaa.update_units(md_calc.unit_conf)
    
    md_calc.strucC = ethane
    
    md_calc.paramC = oplsaa
    
    md_calc.set_ffparam()

.. code:: python

    md_calc.write_data()

