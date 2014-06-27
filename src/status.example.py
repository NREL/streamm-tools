######################################################################
# Load modules
#
# os, sys and subprocess are not picked up from runjobs and must
# be loaded here.
######################################################################

try:
    import os, sys
except:
    print "Error: os and or sys not found"
    print "Check if python is on this system"
    sys.exit(3)

try:
    import subprocess, shlex
except:
    print "Error: subprocess,shlex module not found"
    print "Check 'PYTHONPATH' and/or source bilder *.sh"
    sys.exit(3)

try:
    import runjobs
except:
    print "Error: status module can not find runjobs module"
    print "Check setup.sh path script"
    sys.exit(3)


#
# Derived class overriding default base class method
# This is particular to the application being tested 
#
class ChargeTransJobStatus(runjobs.JobStatus):

    def isJobDone(self, jobdir, jobscript):

        # If output file doesnt not exist, exit as not finished
        fullPath = os.path.join(jobdir,"runData-chgr.dat")
        if (not os.path.exists(fullPath)):
            return False

        # Run external bash line-count
        finishCmd="wc -l " + fullPath
        proc=subprocess.Popen(finishCmd, shell=True,
                              stderr=subprocess.PIPE,
                              stdout=subprocess.PIPE)

        # Success code from cmd above
        return_code = proc.wait()

        # Parse out #-of lines and compare
        for line in proc.stdout:
            stripStr = line.lstrip(' ')
            splstr = stripStr.split(' ')
            fileLines = int(splstr[0])

        # Specific success condition
        if (fileLines == 26):
            return True
        else:
            return False
