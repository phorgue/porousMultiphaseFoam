# -*- coding: utf-8 -*-
## 
## Script for running all validation cases
## \author Pierre Horgue (inspired from MARINE source code and K. Larnier)

# import
from __future__ import with_statement
import os, subprocess, sys, glob, time

# import list_cases
from tutorialsList import tutorials

class testCase:

    #=============================================================================
    # ROUTINE __init__
    #=============================================================================
    def __init__(self, solver, case):
        self.solver = solver
        self.case = case
        self.testDir = f"{solver}-tutorials/{case}"

    #=============================================================================
    # ROUTINE run
    #=============================================================================
    def run(self, preferred_mode):
        print("")
        print(f"Test : {self.solver} {self.case}")
        print("")

        refDir = os.getcwd()
        os.chdir(self.testDir)

        # Check for the existence of run and runParallel scripts
        has_parallel = os.path.exists("./runParallel")
        has_simple = os.path.exists("./run")

        if preferred_mode == "p" and has_parallel:
            method = "parallel"
            script = "./runParallel"
        elif preferred_mode == "s" and has_simple:
            method = "simple"
            script = "./run"
        elif has_parallel:
            print("Preferred method not available. Falling back to parallel.")
            method = "parallel"
            script = "./runParallel"
        elif has_simple:
            print("Preferred method not available. Falling back to simple.")
            method = "simple"
            script = "./run"
        else:
            print("No runnable scripts found for this case.")
            os.chdir(refDir)
            return 1

        print(f"Running in {method} mode...")

        ProcessPipe = subprocess.Popen(
            script, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

        if ProcessPipe:
            stdout, stderr = ProcessPipe.communicate()

        if len(stderr) > 0:
            print("[ ERROR C++ ] " + stderr.decode())

        else:
            foam_exit = "FOAM exiting"
            foam_abort = "FOAM aborting"
            log_files = glob.glob('log.*')
            error_found = False

            for filename in log_files:
                with open(filename, 'r') as foamFile:
                    for line in foamFile:
                        if foam_exit in line or foam_abort in line:
                            error_found = True
                            break
                    if error_found:
                        break

            if error_found:
                print("[ ERROR OpenFOAM ]")
            else:
                print("[ OK ]")

        os.chdir(refDir)
        return 0

#===============================================================================
# PROGRAM Main
#===============================================================================

if __name__ == '__main__':
    print("[ WARNING ] Make sure you have unzipped all the files you need to run your simulations!\n")
    
    parallel = input(str("RUNNING METHOD: \n\tENTER 's' for a simple run\n\tENTER 'p' for a parallel run\nMethod: "))
    
    print("\n========================================================")
    print("                   RUNNING TEST CASES                   ")
    print("========================================================")

    start_time = time.time()
    lap_time = start_time

    for cfg in tutorials:
        solver = cfg["solver"]
        for case_info in cfg["cases"]:
            case = case_info["case"]
            test = testCase(solver, case)
            test.run(parallel)
            print("--- %s seconds ---" % (time.time() - lap_time))
            lap_time = time.time()

    print(" ")
    print("========================================================")
    print("                        FINISHED                        ")
    print("--- %s seconds ---" % (time.time() - start_time))
    print("========================================================")

    sys.exit(0)
