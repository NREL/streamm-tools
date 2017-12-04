**Official docs:** http://streamm.nrel.gov/

The Simulation Toolkit for Renewable Energy and Advanced Materials Modeling (STREAMM) is a python package that generates structures and input files for quantum chemical and molecular dynamics codes.
 STREAMM does not directly conduct simulations rather it is meant to drive quantum chemical and molecular dynamics codes to allow for high-throughput computational analysis of materials.
 The streamm package is written for python2.7 (https://www.python.org/download/releases/2.7/) and incorporates some of the core modules from the pymatgen (http://pymatgen.org/) code.
 

INSTALL
===========


pip install
--------------

To install the streamm package with pip run

    $ pip install streamm

This will create an egg file and place it in your package managment directory of the default version of python on your machine.

Github
--------------

You can also install the streamm package from source by cloning the Github (<https://github.com/>) repo

    $ git clone https://github.com/NREL/streamm-tools
    
By navigating to the package directory
    
    $ cd streamm-tools
    
and running

    $ python setup.py install 

The tests can be run using 

    $ python setup.py test

The package is open-source and can be forked and modified from the repository (github.com/NREL/streamm-tools)


Release Notes
======================

v0.3.4 -- December 2017
----------------------------

* Add P3HT electronic coupling example
* Add export_json/import_json functions to objects
* Cleanup template names 

v0.3.2 -- September 2017
----------------------------

* Add pymatgen (https://github.com/materialsproject/pymatgen) dependency 
* Create tests for each module in the directory tests/
* Split up modules into directories
* Move functions dependent on mpi to util directory

v0.3.0 -- August 2017
----------------------------

* Update the structure of the code to allow `setup.py` installation 


v0.2.0 -- August 28 2015 
----------------------------

* Initial release

NREL (http://www.nrel.gov/) is a National Laboratory of the U.S. Department of Energy,
Office of Energy Efficiency and Renewable Energy, operated by the Alliance for Sustainable Energy, LLC.

Licenses
======================

Streamm license
    
    Copyright 2015 Alliance for Sustainable Energy, LLC
     
    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at
    
        http://www.apache.org/licenses/LICENSE-2.0
    
    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
    
Pymatgen license
    
    The MIT License (MIT)
    Copyright (c) 2011-2012 MIT & LBNL
    
    Permission is hereby granted, free of charge, to any person obtaining a copy of
    this software and associated documentation files (the "Software"), to deal in
    the Software without restriction, including without limitation the rights to
    use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
    the Software, and to permit persons to whom the Software is furnished to do so,
    subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
    FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
    COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
    IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
    CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    
Referencing STREAMM
======================

When referencing the STREAMM toolkit in publications, this website can be cited as.

    Dr. Scott W. Sides, Dr. Travis W. Kemper, Dr. Ross E. Larsen and Dr. Peter Graf.
    "STREAMM (Simulation Toolkit for Renewable Energy and Advanced Materials Modeling),"
    National Renewable Energy Lab, 21 Sept. 2015. <http://github.com/NREL/streamm-tools>.

Also reference the Materials genome project code pymatgen.

    Shyue Ping Ong, William Davidson Richards, Anubhav Jain, Geoffroy Hautier,
    Michael Kocher, Shreyas Cholia, Dan Gunter, Vincent Chevrier, Kristin A.
    Persson, Gerbrand Ceder. *Python Materials Genomics (pymatgen) : A Robust,
    Open-Source Python Library for Materials Analysis.* Computational
    Materials Science, 2013, 68, 314-319. `doi:10.1016/j.commatsci.2012.10.028
    <http://dx.doi.org/10.1016/j.commatsci.2012.10.028>`_
    

Contact us
===========

Email: organicelectronics@nrel.gov
