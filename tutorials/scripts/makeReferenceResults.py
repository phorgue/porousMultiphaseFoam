# -*- coding: utf-8 -*-
## \file checkAndSaveData.py for python 3
## Script for checking tutorials and saving data
## \author Léo Fabrègues & Pierre Horgue

# Import necessary libraries
import shutil, time
from os import path, makedirs

# Import test case definitions
from tutorialsList import tutorials

import re

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
        print(message)
        with open(self.logfile, 'a') as f:
            f.write(message + "\n")  # Append to the log file

#===============================================================================
# PROGRAM Main
#===============================================================================

if __name__ == '__main__':

    # Construct paths
    script_dir = path.dirname(path.abspath(__file__))
    tutorials_dir = path.dirname(script_dir)    
    
    # Extract version information from file
    pmf_dir = path.dirname(tutorials_dir)
    header_file = path.join(pmf_dir, "libraries/general/PMFversion.H")
    with open(header_file, 'r') as file:
        content = file.read()
    match = re.search('const int version_ = (\S+);', content)
    
    if match:
        version = match.group(1)
    else:
        version = "unknown"

    # Initialize logger
    logfile = path.join(tutorials_dir, 'log.plotComparison')
    logfile = "log.references"
    logger = Logger(logfile)

    logger.log("\n========================================================")
    logger.log("      SAVING NEW REFERENCES FROM TUTORIALS       ")
    logger.log("      PMF version is " + version)
    logger.log("========================================================\n")

    # Clean the saveDatas directory if it already exists
    reference_results = path.join(tutorials_dir, "referenceResults")
    if path.exists(reference_results):
        shutil.rmtree(reference_results)
    makedirs(reference_results)
    version_file = open(path.join(reference_results, "version.txt"), 'w')
    version_file.write(version)
    version_file.close()

    # Check and save data for all test cases
    for cfg in tutorials:
        solver = cfg["solver"]
        for id_case in range(len(cfg["cases"])):
            case = cfg["cases"][id_case]["case"]
            result = cfg["cases"][id_case]["result"]
            if not result:
                continue
            source_file = path.join(tutorials_dir, solver+"-tutorials", case, result)
            if path.exists(source_file):
                logger.log("[ OK ] copying file " + source_file)    
                dest_folder= path.dirname(path.join(tutorials_dir, "referenceResults", solver, case, result))
                if not path.exists(dest_folder):
                    makedirs(dest_folder)                    
                logger.log("=====> " + dest_folder)    
                shutil.copy(source_file, dest_folder)

            else:
                logger.log("[ ERROR ] missing file " + source_file)
            
            dest_file = path.dirname(path.join(tutorials_dir, "referenceResults", solver, case, result))

    logger.log("\n========================================================")
    logger.log("                        FINISHED                        ")
    logger.log(f"--- {time.time() - start_time:.2f} seconds ---")
    logger.log("========================================================\n")
