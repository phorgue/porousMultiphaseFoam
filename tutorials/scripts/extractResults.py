# -*- coding: utf-8 -*-
## \file extractResults.py for python 3
## Script for checking tutorials and saving data
## \author Léo Fabrègues (adapted from runTutorials by Pierre Horgue)

# Import necessary libraries
import os, shutil, time

# Import list_cases from the tutorialsList module
from tutorialsList import tutorials as testGroups

start_time = time.time()

class Logger:
    """Class to manage logging in a file while displaying it simultaneously."""
    def __init__(self, logfile):
        self.logfile = logfile
        # Clear the content at the start
        with open(self.logfile, 'w') as f:
            f.write("")                        

    def log(self, message):
        """Method to write and display logs."""
        with open(self.logfile, 'a') as f:
            f.write(message + "\n")  # Append to the log file
        print(message)

class testCase:

    #=============================================================================
    # ROUTINE __init__
    #=============================================================================
    def __init__(self, solver, case, logger):
        self.solver = solver
        self.case = case
        self.testDir = f"{solver}-tutorials/{case}"  # Define the test directory path
        self.logger = logger

    #=============================================================================
    # ROUTINE check_and_save_data
    #=============================================================================
    def check_and_save_data(self):
        self.logger.log(f"\nChecking for postProcessing data in {self.testDir}\n")
        postProcessingDir = os.path.join(self.testDir, 'postProcessing', 'sampleDict')
        
        # If `sampleDict` is missing, check for `probes`
        if not os.path.exists(postProcessingDir):
            postProcessingDir = os.path.join(self.testDir, 'postProcessing', 'probes')
            
        # If `probes` is missing, check for `probesMAP`
        if not os.path.exists(postProcessingDir):
            postProcessingDir = os.path.join(self.testDir, 'postProcessing', 'probesMAP')

        if not os.path.exists(postProcessingDir):
            self.logger.log(f"[ ERROR: missing postProcessing in '{self.testDir}' ]\n")
            return

        # List of subdirectories (which correspond to the time steps)
        time_steps = [d for d in os.listdir(postProcessingDir) if os.path.isdir(os.path.join(postProcessingDir, d))]
        
        # Sort time steps as numbers but keep directory names as strings
        time_steps_sorted = sorted(time_steps, key=lambda x: float(x))
        last_step = time_steps_sorted[-1] if time_steps_sorted else None

        if last_step is None:
            self.logger.log(f"[ ERROR: No valid time step found in '{self.testDir}' ]\n")
            return

        # Destination directory for the saved data
        saveDir = os.path.join(os.getcwd(), "currentResults", f"{self.solver}_{self.case}")
        os.makedirs(saveDir, exist_ok=True)  # Create if it doesn't exist

        sourceDir = os.path.join(postProcessingDir, last_step)
        destDir = os.path.join(saveDir, last_step)

        # Copy files from the last time step to the save directory
        try:
            shutil.copytree(sourceDir, destDir)
            self.logger.log(f"[ OK ] Data from time step {last_step} saved in '{destDir}'\n")
        except Exception as e:
            self.logger.log(f"[ ERROR ] Unable to save data from '{self.testDir}': {e}\n")

# Function to sort saved data by solver
def sort_saved_data_by_solver(currentResultsDir):
    """
    Sorts folders in currentResults by solver, grouping them into subdirectories by solver.
    """
    if not os.path.exists(currentResultsDir):
        print(f"[ ERROR ] The directory {currentResultsDir} does not exist.\n")
        return

    saved_cases = [d for d in os.listdir(currentResultsDir) if os.path.isdir(os.path.join(currentResultsDir, d))]

    for case in saved_cases:
        if "_" not in case:
            print(f"[ ERROR ] Unrecognized case format: {case}\n")
            continue

        solver, caseName = case.split("_", 1)
        solverDir = os.path.join(currentResultsDir, solver)
        os.makedirs(solverDir, exist_ok=True)

        old_path = os.path.join(currentResultsDir, case)
        new_path = os.path.join(solverDir, caseName)

        try:
            shutil.move(old_path, new_path)
            print(f"[ OK ] The case '{case}' has been moved to '{solverDir}' under '{caseName}'\n")
        except Exception as e:
            print(f"[ ERROR ] Unable to move the case '{case}': {e}\n")

#===============================================================================
# PROGRAM Main
#===============================================================================

if __name__ == '__main__':

    logfile = "log.extractResults"
    logger = Logger(logfile)

    logger.log("\n========================================================")
    logger.log("    CHECKING AND EXTRACTING RESULTS FROM TUTORIALS      ")
    logger.log("========================================================\n")

    print("========================================================")
    print("      CHECKING AND EXTRACTING RESULTS FROM TUTORIALS  ")
    print("========================================================\n")

    # Clean the currentResults directory if it already exists
    currentResultsDir = os.path.join(os.getcwd(), 'currentResults')
    if os.path.exists(currentResultsDir):
        shutil.rmtree(currentResultsDir)
    os.makedirs(currentResultsDir)

    # Check and save data for all test cases
    for group in testGroups["tutorials"]:
        solver = group["solver"]
        for case_info in group["cases"]:
            case = case_info["case"]
            test = testCase(solver, case, logger)
            test.check_and_save_data()

    # Sort the data after saving
    sort_saved_data_by_solver(currentResultsDir)
    
    print("\n========================================================")
    print("                        FINISHED                        ")
    print(f"--- {time.time() - start_time:.2f} seconds ---")
    print("========================================================\n")

    logger.log("\n========================================================")
    logger.log("                        FINISHED                        ")
    logger.log(f"--- {time.time() - start_time:.2f} seconds ---")
    logger.log("========================================================\n")
