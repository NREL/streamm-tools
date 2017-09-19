.. _getting_started:

Getting started
***************

The streamm package includes modules to manipulate atomic structures, keep track of force field parameters and write input files for Quantum chemistry and molecular dynamics codes.

First you will create your favorate molecule with a molecular viewer such as `Avogadro <https://avogadro.cc/>`_ and export an .xyz file.
Then read the .xyz file in to a streamm Structure object.::

    import streamm
    fav_mol = Structure('fav')
    fav_mol.read_xyz()
    
    
    
    


