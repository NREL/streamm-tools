tools
=====

Common Tools for all STREAMM projects (private access)

(Note Jan 23 2014: Moved job creation/submission/monitoring python scripts to tools-public)


File structure:

stream/tools/
        src/
                classes
                mpi
                functions
                structures 
                properties
                ( will have Eric set up these files)
        scripts/
                jobmanagment/
                        run_qm
                        run_atomicMD
                structure_generation/
                        molecule_gen.py
                        replicate.py
                force_field/
                        torsion/        
                electron_transfer/

                                
        examples/

        tools/Atom/

Flow:
        All calculations will be tracked via the json file

        0) specify structure
           - molecule_gen.py 
                   - check if exists 
                   - postgresql -> json
                   - json
        2) Run gasphase DFT
           - evaluate conformations 
           - at various levels
           - need to be able to pull structure files form json or form
        database
                - json -> postgresql

        3) Create FF files
           - identify topology
           - fit torisonal potential
        4) Run atomic MD
        5) Analyize MD
           - RDF's
           - Electron transfer
           - json -> postgresql

        6) Create Course grain input
                

************************
Release Goals/Issues:
************************

  Documentation — A writeup of the overall design principles should be included,
                  along with a series of examples and how to run the tests.

  Release details — STREAMM has github, but how exactly does this fit in with the release?
                    Is the github access simply made public? Also, we need to sort out what goes into private repos
                    and public repos. Will the public repos we have setup essentially be like 'tagged' releases?

  Review testing — We should review the current tests and come up with other tests that we need for better coverage

  Portability — How hard are we going to think about this? Its not like we are selling this… but it would be
                good to think about any big barriers to running on different platforms. Again, this comes down to
                actually testing on other platforms… which we should do a minimal amount probably.

  Performance — Serial performance for large systems and timing tests for how the python parallelization is working
