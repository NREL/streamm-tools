.. _installation_instructions:

Installation Instructions
*************************



pip
===

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

    $ python setup.py
    
dependencies will be downloaded and installed and the egg file will placed in your packages directory.
A developer version allowing changes in the source code can be installed using::

    $ python setup.py develop

To access the modules within the package without installing the package sett your `PYTHON_PATH` enviromental variable to include the package source directory `stream/stream`.

The package is open-source and can be forked and modified from the Github repo `<github.com/NREL/streamm>`_.
To merge forked versions of the please contact.




