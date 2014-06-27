#!/usr/bin/env python

"""
Used along with createJobDirs for generating directory structure with job submission
scripts submitting multiple Moab jobs and monitoring queue
"""

try:
    from runjobs import *
except:
    print "Error: runjobs module not found, check PYTHONPATH for runjobs install \n"
    print "Try adding the following to your .bashrc file"
    print "    PYTHONPATH='path-to-runjobs.py':$PYTHONPATH'"
    print "    export PYTHONPATH"
    sys.exit(3)

# Get useful methods
u=Misc()


def createEditors(editType):
    """
    Each element is a file editor object that creates a new
    input file for each new job.
    
    NOTE: this can/should be edited by a user for their own runs
    """

    jobEditorList=list()
    scriptFileEditor=None

    #
    # Creating object that can edit script files
    #
    tagValDictTpl=dict()
    tagValDictTpl['EDIT_rundir'] = "rundirName"
    tagValDictTpl['EDIT_walltime'] = "00:00:30"
    tagValDictTpl['EDIT_nodes'] = 1
    tagValDictTpl['EDIT_ppn']   = 4
    tagValDictTpl['EDIT_procs'] = 4
    scriptFileEditor=FileEditorWithTags(tagValDictTpl)
    tagValDictTpl=None

    #
    # Creating object that can edit input files
    # given a set of rules for param_i ---> param_(i+1)
    # When the getParamSweeps call is used, these are the
    # starting values used for which the 'EDIT_'=... value
    # specified is swept over
    #
    # target distw = 0.005, 0.0125, 0.025, 0.05, 0.10
    #
    tagValDictTpl=dict()
    tagValDictTpl['EDIT_ribsize'] = 20
    tagValDictTpl['EDIT_radius']  = 10.0
    tagValDictTpl['EDIT_distw']   = 0.005
    tagValDictTpl['EDIT_temp' ]   = 273.15
    inFileEditor=FileEditorWithTags(tagValDictTpl)

    params=u.getList(10, 1.0, 26)
    inEditors=inFileEditor.getParamSweeps('EDIT_radius', params)
    jobEditorList.extend(inEditors)

    #
    # Set input file editor explicitly with an
    # assignment string that is parsed to edit to
    # internal dictionary editor 
    # 
    inFileEditor.setDictTplValue("EDIT_distw=0.0125")
    params=u.getList(10, 1.0, 26)
    inEditors=inFileEditor.getParamSweeps('EDIT_radius', params)
    jobEditorList.extend(inEditors)

    #
    # Return one set of editors at a time.
    # This repeats work, but creating these is little work
    # 
    if (editType == 'script'):
        return scriptFileEditor
    elif (editType == 'jobs'):
        return jobEditorList
    else:
        print "Error: createEditors() editType not recognized"
        sys.exit(3)


#
# External run
#
if __name__ == '__main__':

    jobsEditorList=createEditors('jobs')
    scriptFileEditor=createEditors('script')

    print " "
    print "Script editor info \n"
    scriptFileEditor.show()

    sleep(2)
    print " "
    print "Job editors info \n"
    sleep(5)

    for editor in (jobsEditorList):
        editor.show()
