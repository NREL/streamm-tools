.. _templates_example:
  
templates_example
========================
 

In this example, we will write out a NWChem input file based on a
template

.. code:: python

    import streamm

.. code:: python

    import os 
    from pprint import pprint
    from pathlib2 import Path

First, we need to set the location of the templates directory cloned
from http://github.com/NREL/streamm-tools

If you are running this example in the examples directory in the
streamm-tools repo, the TEMPLATE_DIR should look like this

.. code:: python

    EXAMPLE_DIR = os.getcwd()
    
    print(EXAMPLE_DIR)
    
    TEMPLATE_DIR =  os.path.join(EXAMPLE_DIR,'..','templates','')
    
    print(TEMPLATE_DIR)

If not please set the ``TEMPLATE_DIR`` variable to the location of the
templates

We will use the basic nwchem.nw example template

.. code:: python

    temlate_file = 'nwchem.nw'

Create a NWChem calculation object

.. code:: python

    nwchem = streamm.NWChem('ethane_nw_sp')

Read in the ethane structure we creating in the buildingblocks example.

Note::

    If you have not run the buildingblocks_example.ipynb example, please do so to create a `ethane_struc.json` file

.. code:: python

    need_files = ['ethane_struc.json']
    for f in need_files:
        path = Path(f)
        if not path.is_file():
            print("Need to run buildingblocks_example.ipynb")
            os.system("jupyter nbconvert --to python  buildingblocks_example.ipynb")
            os.system("python buildingblocks_example.py")

.. code:: python

    nwchem.strucC.tag = 'ethane'
    nwchem.strucC.import_json()
    print nwchem.strucC.print_properties()

Get the location of the template file

.. code:: python

    template_path =  os.path.join(TEMPLATE_DIR,temlate_file)
    
    print template_path

Read in the template

.. code:: python

    template_line = nwchem.read_lines(template_path)

.. code:: python

    print template_line

Set the properties dictionary to contain the information for our
calculation

.. code:: python

    nwchem.properties['basis'] = '6-31g'
    nwchem.properties['method'] = 'UHF'
    nwchem.properties['charge'] = 0
    nwchem.properties['spin_mult'] = 1
    nwchem.properties['task'] = 'SCF '
    nwchem.properties['coord'] = nwchem.strucC.write_coord()

Do a string replace of the dictionary keys to create an input string

.. code:: python

    input_str = nwchem.replace_keys(template_line,nwchem.properties)
    print input_str

.. code:: python

    file_name = '%s.nw'%(nwchem.tag)
    with open(file_name,"w") as F:
        F.write(input_str)


Easy peasy!
