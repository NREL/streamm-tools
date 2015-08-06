.. _donoracceptorsystems:

.. index:: Gaussian, LAMMPS, DFT, molecular, generator


Molecular Generator: (Donor-Acceptor Systems)
==================================================

Concatenates molecular cply files based on connectivity tags  to generate input
files for `Gaussian <http://www.gaussian.com/>`_ . See :ref:`examples`


Files and directories
-------------------------------------------------------

- BuildingBlocks repos: must contain the directories

   - acceptors   
   - donors            
   - functional_groups spacers           
   - terminals

The backbone of the molecule can be composed of  acceptors, spacers or
donors, and can be  decorated with molecules from the
functional_groups directory and terminated with molecules from the
terminals directory. Run ::

	python ../da_builder/donoracceptorsystems.py -h

to see a list of the options  
