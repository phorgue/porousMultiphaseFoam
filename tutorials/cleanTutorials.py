# -*- coding: utf-8 -*-
## \file cleanTutorials.py
## Script for running all validations cases
## \author Pierre Horgue (inspired from MARINE source code and K. Larnier ;-) )

# import
from __future__ import with_statement
import os, subprocess, sys

# import list_cases
from tutorialsList import tutorials as testCases

class testCase:

    #=============================================================================
    # ROUTINE run
    #=============================================================================
    def __init__(self, solver, case):

        self.solver = solver
        self.case = case
        self.testDir = solver+"-tutorials/"+case

    #=============================================================================
    # ROUTINE run
    #=============================================================================
    def run(self):

        print ""
        print "Test : " + self.solver + " " + self.case
        print ""

        refDir=os.getcwd()

        os.chdir(self.testDir)

        ProcessPipe=subprocess.Popen("./clean", shell=True, \
                                     stdout=subprocess.PIPE,stderr=subprocess.PIPE)

        stdout, stderr = ProcessPipe.communicate()

        print "[ CLEAN ] "

        os.chdir(refDir)

        return 0

#===============================================================================
# PROGRAM Main
#===============================================================================

if __name__ == '__main__':

    print "========================================================"
    print "                   CLEANING TEST CASES                  "
    print "========================================================"

    for case in testCases:
        test = testCase(case["solver"],case["case"])
        test.run()

    print " "
    print "========================================================"
    print "                        FINISHED                        "
    print "========================================================"

    sys.exit(0)
