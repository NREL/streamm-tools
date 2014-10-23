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
.. automodule:: test_particleContainer
   :members:

.. automodule:: test_angleContainer
   :members:



