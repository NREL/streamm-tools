.. _tools_tests:

.. index:: regression tests, unit tests, test scripts


*******************************
Unit tests
*******************************


Running tests
========================

The tests for the STREAMM tools is in a separate Github repository https://github.nrel.gov/streamm/tools-tests hosted by the National
Renewable Energy Lab (NREL). The tools repo can be downloaded using the linux command line ::

     git clone https://github.nrel.gov:streamm/tools-tests.git

If the tools repo has been configured correctly the check.sh script
can be used to run the units tests. The help info for check.sh is
output from ::

    check.sh

To run a particular test ::

    check.sh 'test-name.py'

and to run all tests ::

    check.sh all



Test descriptions
========================

.. automodule:: test_angleContainer
   :members: main

.. automodule:: test_bondContainer
   :members: main

.. automodule:: test_checkTypes
   :members: main

.. automodule:: test_dihedralContainer
   :members: main

.. automodule:: test_dumpLammps
   :members: main

.. automodule:: test_elements
   :members: main

.. automodule:: test_ffparameters
   :members: main

.. automodule:: test_n2_mpiBase
   :members: main

.. automodule:: test_nX_mpiBase_splitData
   :members: main

.. automodule:: test_particleConstructors
   :members: main

.. automodule:: test_particleContainer
   :members: main

.. automodule:: test_particleSetInfo
   :members: main

.. automodule:: test_ptclContainerConstructor
   :members: main

.. automodule:: test_searchTags
   :members: main

.. automodule:: test_simulation
   :members: main

.. automodule:: test_strucAdd
   :members: main

.. automodule:: test_strucAdd2
   :members: main

.. automodule:: test_strucAddBig
   :members: main

.. automodule:: test_strucCompressID
   :members: main

.. automodule:: test_strucCompressIDWithAngles
   :members: main

.. automodule:: test_strucDumpSave
   :members: main

.. automodule:: test_strucEmpty
   :members: main

.. automodule:: test_strucEmpty2
   :members: main

.. automodule:: test_strucSetPtclPos
   :members: main

.. automodule:: test_strucWithAngles
   :members: main

.. automodule:: test_strucWithDihedrals
   :members: main

.. automodule:: test_subStructure
   :members: main

.. automodule:: test_subStructureWithAngles
   :members: main

.. automodule:: test_subStructureWithDihedrals
   :members: main
