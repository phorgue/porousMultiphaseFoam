# -*- coding: utf-8 -*-
## \file cleanTutorials.py for python 3
## Script for cleaning all validation cases
## \author Pierre Horgue (inspired from MARINE source code and K. Larnier)

# import
from __future__ import with_statement
import os, subprocess, sys, glob, shutil

# Import test case definitions
from tutorialsList import tutorials

class testCase:

    #=============================================================================
    # ROUTINE __init__
    #=============================================================================
    def __init__(self, tutorials_dir, solver, case):
        self.solver = solver
        self.case = case
        self.testDir = os.path.join(tutorials_dir, f"{solver}-tutorials/{case}")

    #=============================================================================
    # ROUTINE run
    #=============================================================================
    def run(self):
        print("")
        print(f"Cleaning Test : {self.solver} {self.case}")
        print("")

        os.chdir(self.testDir)

        ProcessPipe = subprocess.Popen(
            "./clean", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        stdout, stderr = ProcessPipe.communicate()

        print("[ CLEAN ] ")

        return 0

#===============================================================================
# PROGRAM Main
#===============================================================================

if __name__ == '__main__':

    # Construct paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    tutorials_dir = os.path.dirname(script_dir)

    print("========================================================")
    print("                   CLEANING TEST CASES                  ")
    print("========================================================")

    for cfg in tutorials:
        solver = cfg["solver"]
        for case_info in cfg["cases"]:
            case = case_info["case"]
            test = testCase(tutorials_dir, solver, case)
            test.run()

    # Create directory for figures and clean if already exists
    figures_dir = os.path.join(tutorials_dir, "validationFigures")
    if os.path.exists(figures_dir):
        print("\nRemove validationFigures dir")
        shutil.rmtree(figures_dir)

    os.chdir(tutorials_dir)
    for f in glob.glob("log.*"):
        os.remove(f) 

    print(" ")
    print("========================================================")
    print("                        FINISHED                        ")
    print("========================================================")

    sys.exit(0)
