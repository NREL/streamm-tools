#!/usr/bin/env python

######################################################################
#
# S. Sides  (NREL)
# 
# Script for submitting multiple Moab jobs and monitoring queue
#
######################################################################
try:
    from runjobs import *
except:
    print "Error: runjobs module not found, check PYTHONPATH for runjobs install"
    print "Try adding the following to your .bashrc file"
    print "   PYTHONPATH='path-to-runjobs.py':$PYTHONPATH'"
    print "   export PYTHONPATH"
    sys.exit(3)

try:
    import numpy
except:
    print "Error: numpy module not found, check PYTHONPATH for numpy install or module version"
    print "Try adding the following to your .bashrc file"
    print "   PYTHONPATH='path-to-runjobs.py':$PYTHONPATH'"
    print "   export PYTHONPATH"
    sys.exit(3)

# Get useful methods
u=Misc()
#############################################################################


#############################################################################
# Command line option parse

parser = OptionParser()
usage = """

%prog [options]

     Manages multiple job submission, queue monitoring and
     checks for completed simulations

     ##############
     Status files:
     ##############
          jobs-submitJobs.log --> status info from script
       jobs-submit-failed.log --> failed job submit info (location, scriptname)/line
      jobs-submit-success.log --> successful job submit info (location, scriptname)/line
            jobs-run-done.log --> job code finished list
         jobs-run-notdone.log --> job code not finished list

     NOTE: This script can submit MANY jobs. If USERNAME needs to kill
           these jobs before they run use 'mjobctl -c -w user=USERNAME'

     NOTE: The information in these files is always from the last submit command.
           It will not track the history of submit results.


     ##########################
     Simulation finish checks:
     ##########################
        If user provides a specific "code finished" condition, when the script is re-run
        only jobs that are not finished will be resubmitted (job is restarted for initial run)

        The method to be edited for different "code finished" conditions is
         
           def isJobDone(jobdir, jobscript):
               ...
               ...
               return (True/False)

        where 'jobdir' is the path to the run directory and jobscript is the name
        of the job submission script. the method needs to return
        'True' if the code finished successfully and 'False' otherwise
        The derived class from from the JobStatus class overrides the isJobDone method
        with application specific code.
        status.py should be in PYTHONPATH and the correct object should be chosen below.
        Look for JobStatus: note in submitJobs.py script

"""
parser.set_usage(usage)

parser.add_option("-s","--search",
                  dest="searchPattern",
                  default='.mb',
                  type="string",
                  help="Search pattern for submission scripts (eg .slurm, .mb) \n")

parser.add_option("-t","--topdir",
                  dest="topdir",
                  type="string",
                  help="Root directory to begin script search \n")

parser.add_option("-e","--editstring",
                  dest="editstring",
                  default=None,
                  type="string",
                  help="Edit found script files eg. -e 'search-string --> target-string' \n")

parser.add_option("-d","--dryrun",
                  dest="dryrun",
                  action='store_true',
                  default=False,
                  help="Runs through script without submitting \n")

parser.add_option("-c","--checkdone",
                  dest="checkdone",
                  action='store_true',
                  default=False,
                  help="Check number of jobs completed [need to set -t option] \n")

parser.add_option("-m","--maxAttempts",
                  dest="maxAttempts",
                  type="int",
                  default=1,
                  help="Maximum #-attempts to submit job (default=1) \n")

# Acessing options
(options,args) = parser.parse_args()
searchPattern  = options.searchPattern
topdir         = options.topdir
dryrun         = options.dryrun
checkdone      = options.checkdone
maxAttempts    = options.maxAttempts
editstring     = options.editstring

# Check if options set
if (topdir==None):
    options_set=False
elif (checkdone):
    options_set=False  # Setting false so main section skipped
    print " "
    print "Will check for completed jobs in ", topdir, " \n"
elif (editstring):
    options_set=False
    print " "
    print "Will search script files for edit string and replace \n"
else:
    # Check top direcotry which must be present
    options_set=True
    print "Top directory ", topdir
    if ( not os.path.exists(topdir) ):
        print "Top directory (specify with -t) doesn't exist \n"
        sys.exit(1)
#############################################################################




####################################################################
#
# Runs recursive glob on directory structure starting at
#   'td'  -- top directory to search for scripts
#   'ptn' -- script name search string
#
# Checks paths and formats paths, file names
# Runs 'finish' condition check and returns in tuple
#
####################################################################
def getJobInfo(td, ptn):

    global stdIO
    global u

    # Local lists
    jobsInfo=[]
    fullPathNames=[]

    # Check paths
    tdExists = os.path.exists(td)
    if (not tdExists):
        stdIO.prt("Directory " + td + " does not exist \n")
        sys.exit(1)

    # Main glob call, adds wildcard
    wldPtn = "*" + ptn
    files=u.recursive_glob(td, wldPtn)

    # Parse paths and script names
    # eg f = test-runs/runDir1/runDir-1.1//job.mb
    for f in files:
        job=[]
        jobPath = os.path.dirname(f)    # test-runs/runDir1/runDir-1.1
        jobScript = os.path.basename(f) # job.mb
        # Determine is job finished
#        jobDone = isJobDone(jobPath, jobScript)
        jobDone = js.isJobDone(jobPath, jobScript)

        job.append(jobPath)             # .
        job.append(jobScript)           # ..
        job.append(jobDone)             # job = (test-runs/runDir1/runDir-1.1, job.mb, True/False)
        jobsInfo.append(job)            # jobsInfo(job,....)

    # Check for found jobs and exit if empty
    if ( len(jobsInfo) == 0 ):
        stdIO.prt("No job submission scripts found with " + ptn)
        stdIO.prt("Check the -s option and spelling \n")
        parser.print_help()
        sys.exit(1)

    return jobsInfo

##########################################################
#
# Submit job to queue and record stdout and stderr
# info. This method is fairly general, but only tested
# for the msub command. Therefore, the msub command is
# explictly constructed here. This should be executed
# once inside job submission directory
#
# scriptName -- name of job submission script (eg job.sh)
#
##########################################################
def manageJob(scriptName):

    global stdIO
    global failJobIO
    global goodJobIO

    # Set submit command
    msubcommand = "qsub " + scriptName
    
    # If dry run test return with 'fail' code
    if dryrun:
        stdIO.prtFile("Would execute: " + msubcommand)
        return 1
    else:
        stdIO.prtFile("Executing: " + msubcommand)

    # Runs submission command and grabs output
    proc = subprocess.Popen(msubcommand, shell=True,
                            stderr=subprocess.PIPE,
                            stdout=subprocess.PIPE)

    # Success code from msub command
    return_code = proc.wait()

    # Current max queue submit rate (just for safety)
    sleep(0.20)

    # Read from pipes
    # Note: Stripping out blank output from moab
    for line in proc.stdout:
        if (len(line.rstrip()) != 0):
            goodJobIO.prtFile(line.rstrip())
            stdIO.prtFile("jobID: " + line.rstrip())

    for line in proc.stderr:
        if (len(line.rstrip()) != 0):
            thisdir = os.getcwd()
            failJobIO.prtFile(thisdir + " " + line.rstrip())
            failJobIO.flush()
            stdIO.prtFile("Failed w/ " + line.rstrip())

    return return_code


#####################################################################
#                              main                                 #
#####################################################################

# Globals
topDebug=False     # Developer debug codes
queueJobLimit=1000 # HPC queue limit
stdIO=None
failJobIO=None
goodJobIO=None
goodJobs=0
failJobs=0

#
# JobStatus:
# Pick the correct derived class supplying the
# isJobDone method depending on code application
#
js=JobStatus()
try:
    from status import GaussianJobStatus
    js=GaussianJobStatus()
    print "Gaussian job status module found"
except:
    print "status for Gaussian not found"
try:
    from status import ChargeTransJobStatus
    js=ChargeTransJobStatus()
    print "ChargeTrans job status module found"
except:
    print "status for ChargeTrans not found"


# Main submit jobs section
if (options_set):

    # List of 2-tuple (job location, job script)
    jobs = getJobInfo(topdir, searchPattern)
    numJobs = len(jobs)

    # Create IO objects, class cleans up on exit
    stdFname = os.path.join(topdir,"jobs-submitJobs.log")
    stdIO=IO(stdFname, topDebug)
    stdIO.prtFile("submitJobs.py Log \n")

    failFname = os.path.join(topdir,"jobs-submit-failed.log")
    failJobIO=IO(failFname, topDebug)
    failJobIO.prtFile("Failed job submission info \n")

    goodFname = os.path.join(topdir,"jobs-submit-success.log")
    goodJobIO=IO(goodFname, topDebug)
    goodJobIO.prtFile("Successful submission jobIDs \n")

    doneFname = os.path.join(topdir,"jobs-run-done.log")
    doneJobIO=IO(doneFname, topDebug)
    doneJobIO.prtFile("Jobs run finished \n")

    notFname = os.path.join(topdir,"jobs-run-notdone.log")
    notJobIO=IO(notFname, topDebug)
    notJobIO.prtFile("Jobs run failed/not-complete \n")

    # Extra status info
    stdIO.prt("Maximum submit attempts/job = " + str(maxAttempts))

    # Number of scripts found
    stdIO.prt("Found " + str(numJobs) + " jobs \n")

    for j in jobs:

        jdir  = j[0]             # Job path
        jscr  = j[1]             # Job script
        jdone = j[2]             # Job status

        # If job is complete move on
        if (jdone):
            doneJobIO.prt(jdir + " is complete... going to next job")
            continue
        else:
            notJobIO.prt(jdir + " is not complete, will (re)/submit")

        currdir = os.getcwd()    # 'this' script location
        os.chdir(jdir)           # Change to job script location

        result = manageJob(jscr) # First attempt
        submitAttempts=1         # Set job state parameters

        # Job-resubmission logic
        while ( (result != 0) and (submitAttempts < maxAttempts)):
            print "Job submission ", submitAttempts, " failed, re-attempting ..."
            submitAttempts = submitAttempts + 1
            result = manageJob(jscr)

        # Job-resubmission logic
        if ( (result !=0) and (submitAttempts == maxAttempts)):
            print "Exceeded max attempt limit of ", maxAttempts, " ...moving on"
        
        os.chdir(currdir)        # Change back to original directory

        # Print/record results of job submission
        if (result==0):
            goodJobs=goodJobs+1
            stdIO.prt(jdir + "/" + jscr + " successfully in queue \n")
            stdIO.flush()
            # HPC imposed limit
            if goodJobs > queueJobLimit:
                stdIO.prt("HPC queue job limit " + str(queueJobLimit) + " reached.")
                break
        else:
            failJobs=failJobs+1
            stdIO.prt(jdir + "/" + jscr + " failed to get in queue \n")
            stdIO.flush()

    # Report summary
    stdIO.prt(" ")
    stdIO.prt("Successful job submissions --> " + str(goodJobs))
    stdIO.prt("    Failed job submissions --> " + str(failJobs))
    stdIO.prt(" ")



# Only check which jobs are done
elif (checkdone):

    completedJobs = 0
    unfinishedJobs = 0

    # List of 2-tuple (job location, job script)
    jobs = getJobInfo(topdir, searchPattern)
    numJobs = len(jobs)

    for j in jobs:

        jdir  = j[0]             # Job path
        jscr  = j[1]             # Job script
        jdone = j[2]             # Job status

        # If job is complete move on
        if (jdone):
            completedJobs = completedJobs + 1
        else:
            unfinishedJobs = unfinishedJobs + 1

    print "Jobs: "
    print "   completed  = ", completedJobs
    print "   unfinished = ", unfinishedJobs
    print "   total      = ", numJobs
    print " "



# Use editstring to change all submit script (for matching line)
elif (editstring):

    # Parse -e option line string
    twostrings = string.split(editstring, "-->")
    if (len(twostrings) != 2):
        print "Error with -e editstring format"
        print "Make sure '-->' delimiter in place \n"
        sys.exit(1)

    # Further parse -e option for search/replace strings
    searchStr  = twostrings[0].strip(' ')
    replaceStr = twostrings[1].strip(' ')
    print "Editing all found script files"
    print "  search-string  = ", searchStr
    print "  replace-string = ", replaceStr
    print " "

    # List of 2-tuple (job location, job script)
    jobs = getJobInfo(topdir, searchPattern)
    numJobs = len(jobs)

    for j in jobs:

        jdir  = j[0]             # Job path
        jscr  = j[1]             # Job script
        jdone = j[2]             # Job status

        # Construct full path to script file
        scriptPath = jdir + "/" + jscr
        print "Found ",scriptPath

        # If dryrun just print info and dont edit
        if (dryrun):
            print "... would edit ", scriptPath
            continue

        # Do simple 'grep' for search string
        outFile="tmp.out"
        fobj=open(scriptPath, 'r+')
        gobj=open(outFile,    'w')
        for line in fobj:
            if searchStr in line:
                print "... editing ", scriptPath
                gobj.write(replaceStr)
                gobj.write("\n")
            else:
                gobj.write(line)

        # Close local files and
        # copy over tmp edited file to script file
        fobj.close()
        gobj.close()
        shutil.copy2(outFile, scriptPath)

    # Cleanup
    os.remove(outFile)


# If no options, print help message
else:
    parser.print_help()
