# -*- coding: utf-8 -*-
## \file runTutorials.py for python 3
## Script for running all validations cases
## \author Pierre Horgue (inspired from MARINE source code and K. Larnier)

# import
from __future__ import with_statement
import os, subprocess, sys, glob, time

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

        print("")
        print("Test : " + self.solver + " " + self.case)
        print("")

        refDir=os.getcwd()

        os.chdir(self.testDir)

        ProcessPipe=subprocess.Popen("./run", shell=True, \
                                     stdout=subprocess.PIPE,stderr=subprocess.PIPE)

        stdout, stderr = ProcessPipe.communicate()

        if len(stderr) > 0:

            print("[ ERROR C++ ] " + stderr)

        else:

            foam_exit = "FOAM exiting"
            foam_abort = "FOAM aborting"
            log_files = glob.glob('log.*')
            for filename in log_files:
                foamFile = open(filename, 'r')
                error_found = False
                for lines in foamFile:
                    if foam_exit in lines:
                        error_found = True
                        break
                    if foam_abort in lines:
                        error_found = True
                        break
                    
            if error_found:
                print("[ ERROR OpenFOAM ] ")
            else:
                print("[ OK ] in ")

            os.chdir(refDir)

        return 0

#===============================================================================
# PROGRAM Main
#===============================================================================

if __name__ == '__main__':

    print("========================================================")
    print("                   RUNNING TEST CASES                   ")
    print("========================================================")

    start_time = time.time()
    lap_time = start_time
    for case in testCases:
        test = testCase(case["solver"],case["case"])
        test.run()
        print("--- %s seconds ---" % (time.time() - lap_time))
        lap_time = time.time()

    print(" ")
    print("========================================================")
    print("                        FINISHED                        ")
    print("--- %s seconds ---" % (time.time() - start_time))
    print("========================================================")

    sys.exit(0)
