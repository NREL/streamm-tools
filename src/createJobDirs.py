#!/usr/bin/env python

"""
For generating directory structure with job submission
scripts submitting multiple Moab jobs and monitoring queue
"""

try:
    from runjobs import *
except:
    print "Error: runjobs module not found, check PYTHONPATH for runjobs install \n"
    print "Try adding the following to your .bashrc file"
    print "    PYTHONPATH='path-to-runjobs.py':$PYTHONPATH'"
    print "    export PYTHONPATH"
    import sys
    sys.exit(3)

# Get useful methods
u=Misc()
######################################################################


#############################################################################
# Command line option parse

parser = OptionParser()
usage = """

%prog [option]

     Manages multiple job creation w/parameter scans
     suitable for submitting through a queueing system
     
     User must set up a 'repo' directory with:
        * an executable
        * an input file 'template'
        * a PBS/Moab queue script 'template'

     Assumptions: each job to be submitted corresponds to
     a separate directory with an edited input file and if
     necessary an edited script file. The script file uses
     the executable and edited input file template in the
     appropriate way for the specific code.

     NOTE: the term 'template' refers to a file that can
     be edited (or otherwise transformed) by the runjobs
     python module.

     NOTE: An editors.py file must be in PYTHONPATH. This
     file contains methods to create the editor objects for
     scripts and input files. This separate file allows
     same create jobs and submit jobs .py scripts for
     different input parameters.
"""
parser.set_usage(usage)

parser.add_option("-t","--topRunDir",
                  dest="topRunDir",
                  type="string",
                  default=os.getcwd(),
                  help="Root directory to begin job creation \n")

parser.add_option("-v","--verbosity",
                  dest="jverbosity",
                  action='store_true',
                  default=False,
                  help="Run with verbose info \n")

parser.add_option("-d","--dryrun",
                  dest="dryrun",
                  action='store_true',
                  default=False,
                  help="Runs through script parameters \n")

# Acessing options
numargs=len(sys.argv)-1
(options,args) = parser.parse_args()
topRunDir  = options.topRunDir
jverbosity = options.jverbosity
dryrun     = options.dryrun
#############################################################################


def showInfo():
    """
    Print parameter settings to screen
    """

    print " "
    print "************************* "
    print "Job creation settings:    "
    print "************************* "
    print " "
    print "  To change these, edit this script"
    print "    # Path to repo for common files             - repoDirPath = ", repoDirPath
    print "    # Name of executable binary in repo         - execName    = ", execName
    print "    # Name of input file (template) in repo     - inputName   = ", inputName
    print "    # Name of script for PBS (template) in repo - scriptName  = ", scriptName
    print "    # Auxillary files to copy                   - auxFileList = ", auxFileList
    print "   "
    print "  To change these defaults below, set options (see help, run with -h)"
    print "    # Top-level dir where runs will be created  - topRunDir   = ", topRunDir
    print "    # Verbosity switch                          - jverbosity  = ", jverbosity
    print " "
    import editors
    print "Compiled job edit module found at ", editors.__file__
    print " "


def runCreateJobsDirs(jobj, jobEditorList, scriptEditor):
    """
    Job directory names default to "run-00001, run-00002..."

    Args:
        jobj          -- is a properly initialized Job() object
        jobEditorList -- is a list of input file edit objects
        scriptEditor  -- file editor object for the submit script

    Returns: None
    """
    
    global u

    # Set basic script editor for all jobs
    jobj.setScriptEditor(scriptEditor)

    #
    # Loop over input file editor objects created above
    # Each input file object represents a new job because the editor
    # will be used to create a new input file in each new run directory
    #
    for ijob, jobInputEditor in enumerate(jobEditorList):
    
        # Construct next job setup object from template
        j=copy.copy(jobj)
    
        # Default run directory name
        jobNumStr = u.getTagForSorting(ijob+1)
        runDir="run-" + jobNumStr
        fullRunDir = os.path.join(topRunDir, runDir)

        # Set script edit so run directory is inserted
        scriptEditor.setDictValue("EDIT_rundir="+runDir)

        j.setRunDir(fullRunDir)             # Set run directory pathname
        j.setInputEditor(jobInputEditor)    # Set job object input  file editor
        j.createJobRun()                    # Create run space with file from repo
        j.changeInputFile()                 # Use file editor object on inputfile
        j.changeScriptFile()                # Use file editor object on scriptfile





##################################################################
#
# Set host specific paths /commands
#
##################################################################
hostname=os.getenv('HOSTNAME')

if ( "login" in hostname or "n0" in hostname):
    prospect_exe="prospect-intel"
elif ( "stc" in hostname ):
    prospect_exe="prospect-gcc"
else:
    print "Hostname not recognized"
    sys.exit(0)
##################################################################



########################################################################################
#
# USER-Editable: Main parameters
# 
########################################################################################
repoDirPath='./jobrun-template'      # Path to repo for common files
# execName='prospect-intel'          # Name of executable binary in repo
execName=prospect_exe                # Name of executable binary in repo
inputName='prospect.in'              # Name of input file (template) in repo
scriptName='script-chg.pbs'          # Name of script for PBS (template) in repo
auxFileList=[
'raw2vtk.py',
'prospect-viz.in'
]                                    # Names of auxillary files to copy to run dirs
########################################################################################


####################################################################################################
#
# Run method to create input file
# editors and script file editor
#
####################################################################################################

try:
    from editors import *
except:
    print "Error: editors module not found, check PYTHONPATH for editors install \n"
    print "Try adding the following to your .bashrc file"
    print "    PYTHONPATH='path-to-editors.py':$PYTHONPATH'"
    print "    export PYTHONPATH"
    sys.exit(3)

__name__ = 'internal'
jobsEditorList=createEditors('jobs')
scriptFileEditor=createEditors('script')
####################################################################################################


#
# Just print help message
#
if (numargs==0):
    parser.print_help()

#
# Give dry run info
#
elif (dryrun):

    showInfo()
    if (jverbosity):
        print " "
        print "Will output edit info \n"
        sleep(2)
        for jobEditor in (jobsEditorList):
            jobEditor.show()

#
# Create job directories ('main')
#
else:

    j0=Job(repoDirPath, execName, inputName, scriptName, auxFileList, verbose=jverbosity)  # Creating initial job object
    runCreateJobsDirs(j0, jobsEditorList, scriptFileEditor)                                # Main (changing script editor in method)

    # Copy the editor.py module used to create jobs
    # in top run directory
    if (not os.path.exists("editors.py")):
        print " "
        print "Editors module not found"
        print sys.exit(1)
    else:
        shutil.copy2("editors.py", topRunDir)
        print " "
        print "Editors module copied to ", topRunDir

    print " "
    print "Job directories w/scripts created. Submit jobs by:"
    print "  1. msub each script 'by hand'"
    print "  2. use submitJobs.py script to automate"
    print " "

# Cleaning auxillary files
os.system('rm -rf ./*.pyc')

