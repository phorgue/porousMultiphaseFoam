# -*- coding: utf-8 -*-
## \file runTutorials.py
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

        ProcessPipe=subprocess.Popen("./run", shell=True, \
                                     stdout=subprocess.PIPE,stderr=subprocess.PIPE)

        stdout, stderr = ProcessPipe.communicate()

        if len(stderr) > 0:

            print "[ ERROR C++ ] " + stderr

        else:

            find_str = "FOAM exiting"
            foamFile = "log."+self.solver
            with open(foamFile, "r") as f:
                f.seek(0, 2)
                fsize = f.tell()
                f.seek (max (fsize-1024, 0), 0)
                lines = f.readlines()

            lines = lines[-3:]

            if find_str + "\n" in lines:

                print "[ ERROR OpenFOAM ] "

            else:

                print "[ OK ]"

            os.chdir(refDir)

        return 0

#===============================================================================
# PROGRAM Main
#===============================================================================

if __name__ == '__main__':

    print "========================================================"
    print "                   RUNNING TEST CASES                   "
    print "========================================================"

    for case in testCases:
        test = testCase(case["solver"],case["case"])
        test.run()

    print " "
    print "========================================================"
    print "                        FINISHED                        "
    print "========================================================"

    sys.exit(0)
