.. _installation_instructions:

Installation Instructions
*************************

The streamm package is writen for python2.7 and is dependent on some core modules from the pymatgen code. 

pip intall
==========

The pip command is the standard method for installing python packages, and is included in most distributions of python. 
To install the streamm package with pip run::

    $ pip install streamm

This will create an egg file and place it in your package managment directory of the default version of python on your machine.
The pip command can be installed on MacOS/PC/Linux according to `<https://pip.pypa.io/en/stable/installing/>`_.

Github
======

The streamm package can be also installed from source by cloning the github repo::

    $ git clone https://github.com/NREL/streamm
    
By navigating to the package directory::
    
    $ cd streamm
    
and running::

    $ python setup.py install 
    
dependencies will be downloaded and installed and the egg file will placed in your packages directory.
The tests can be ran using ::

    $ python setup.py test

The package is open-source and can be forked and modified from the Github repo `<github.com/NREL/streamm>`_.