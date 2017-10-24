.. _templates:
  
templates
===============
 

In this example, we will write out a NWChem input file based on a
template

.. code:: python

    import streamm

.. code:: python

    from __future__ import division, unicode_literals

.. code:: python

    import os 
    from pprint import pprint

First, we need to set the location of the templates directory cloned
from http://github.com/NREL/streamm-tools

If you are running this example in the examples directory in the
streamm-tools repo, the TEMPLATE\_DIR should look like this

.. code:: python

    EXAMPLE_DIR = os.getcwd()
    
    print(EXAMPLE_DIR)
    
    TEMPLATE_DIR =  os.path.join(EXAMPLE_DIR,'..','templates','')
    
    print(TEMPLATE_DIR)


.. parsed-literal::

    /Users/rlarsen/Development/streamm-tools/examples
    /Users/rlarsen/Development/streamm-tools/examples/../templates/


If not please set the ``TEMPLATE_DIR`` variable to the location of the
templates

We will use the basic nwchem.nw example template

.. code:: python

    temlate_file = 'nwchem.nw'

Create a NWChem calculation object

.. code:: python

    nwchem = streamm.NWChem('ethane_nw_sp')

Read in the ethane structure we creating in the buildingblocks example.

.. Note::

::

    If you have not run the buildingblocks.ipynb example, please do so the create a `.xyz` file

.. code:: python

    nwchem.strucC.tag = 'ethane'
    nwchem.strucC.read_xyz()
    print nwchem.strucC.n_particles


.. parsed-literal::

    8


Get the location of the template file

.. code:: python

    template_path =  os.path.join(TEMPLATE_DIR,temlate_file)
    
    print template_path


.. parsed-literal::

    /Users/rlarsen/Development/streamm-tools/examples/../templates/nwchem.nw


Read in the template

.. code:: python

    template_line = nwchem.read_lines(template_path)


::


    ---------------------------------------------------------------------------

    AttributeError                            Traceback (most recent call last)

    <ipython-input-15-f054b2c3cca3> in <module>()
    ----> 1 template_line = nwchem.read_lines(template_path)
    

    AttributeError: 'NWChem' object has no attribute 'read_lines'


.. code:: python

    print template_line


.. parsed-literal::

    start test
     geometry GEOM units angstroms NOCENTER NOAUTOZ NOAUTOSYM
    <coord>end
    
    
     BASIS 
     * LIBRARY <basis>
     end 
     SET geometry  GEOM 
     CHARGE  <charge>
     SCF 
     NOPEN 0
     <method> 
     SINGLET
     maxiter 100
     end 
     TASK <task>
    
    


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


.. parsed-literal::

    start test
     geometry GEOM units angstroms NOCENTER NOAUTOZ NOAUTOSYM
         C       1.34000000      -0.00000000       0.00000000 
         H       1.74000000      -0.00000000      -1.13137084 
         H       1.74000000       0.97979589       0.56568542 
         H       1.74000000      -0.97979589       0.56568542 
         C       0.00000000       0.00000000       0.00000000 
         H      -0.40000000       0.00000000       1.13137084 
         H      -0.40000000      -0.97979589      -0.56568542 
         H      -0.40000000       0.97979589      -0.56568542 
    end
    
    
     BASIS 
     * LIBRARY 6-31g
     end 
     SET geometry  GEOM 
     CHARGE  0
     SCF 
     NOPEN 0
     UHF 
     SINGLET
     maxiter 100
     end 
     TASK SCF 
    
    


.. code:: python

    file_name = '%s.nw'%(nwchem.tag)
    with open(file_name,"w") as F:
        F.write(input_str)


Easy peasy!
