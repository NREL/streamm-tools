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
                


