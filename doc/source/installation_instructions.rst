.. _installation_instructions:

Installation Instructions
********************************

`Do we want to specify Mac/PC/Linux?`

pip
===

The pip command is the standard method for installing python packages, and is included in most distributions of python. 
To install the streamm package with pip run::

    $ pip install streamm

This will create an egg file and place it in your package managment directory of the default version of python on your machine.
The pip command can be installed on MacOS/PC/Linux according to `<https://pip.pypa.io/en/stable/installing/>`.

Github
======

The streamm package can be also installed from source by cloning the github repo::

    $ git clone https://github.com/NREL/streamm
    
By navigating to the package directory::
    
    $ cd streamm
    
and running::

    $ python setup.py
    
dependencies will be downloaded and installed and the egg file will placed in your packages directory.
The modules within the package directory can be accesses by setting your `PYTHON_PATH` enviromental variable to include `stream/stream`.


