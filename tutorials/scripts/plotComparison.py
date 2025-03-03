#!/usr/bin/env python3
# -*- coding: utf-8 -*-
## Script for plotting current results VS reference results
## \author Léo Fabrègues & Pierre Horgue

# Imports
from __future__ import with_statement
import glob, shutil, csv
from os import path, makedirs
import numpy as np
import pylab as plt
plt.rc('text', usetex=False)
import re

# Import test case definitions
from tutorialsList import tutorials

# Logger class to handle logging
class Logger:
    """Class to manage logging in a file while displaying it simultaneously."""
    def __init__(self, logfile):
        self.logfile = logfile

    def log(self, message):
        """Write message to the log file and print it."""
        print(message)
        with open(self.logfile, 'a') as f:
            f.write(message + "\n")

# Function to determine the PMF version
def get_pmf(currentDirectory):
    # Extracts the directory path for PMF version
    pmfVersionDirectory = currentDirectory.split('/')[1:]
    res = []
    
    for dir in pmfVersionDirectory:
        if dir == 'porousMultiphaseFoam':
            res.append(dir)
            break
        else:
            res.append(dir)
    
    pmfVersionDirectory = path.join('/', *res, 'solvers')
    target_file = 'headerPMF.H'
    pmfVersion = 'current_version'

    # Extract version information from file
    def extract_pmf_version(file_path):
        with open(file_path, 'r') as file:
            content = file.read()
        match = re.search(r'Current version of PMF is (\S+)', content)
        return match.group(1) if match else None

    # Verify and extract PMF version
    file_path = path.join(pmfVersionDirectory, target_file)
    if path.exists(file_path):
        pmfVersion = extract_pmf_version(file_path)
        return pmfVersion[:(len(pmfVersion)-1)]
    else:
        print('[ ERROR ] No Header Found')
        
    return pmfVersion


# Function to process .csv files and generate plots
def extract_XY(main_dir, config):

    # Read data from CSV file
    def read_csv_data(file_path):
        data_x = []
        data_y = []
        dataLabels = []

        with open(file_path, 'r') as file:
            csvreader = csv.reader(file)
            for i, line in enumerate(csvreader):
                if i == 0:
                    dataLabels = line[1:]
                else:
                    data_x.append(float(line[0]))
                    data_y.append([float(value) for value in line[1:]])

        return np.array(data_x), np.array(data_y).T, dataLabels


    file_path = path.join(main_dir, case, config["result"])
    resultPresent = path.isfile(file_path)

    if resultPresent:
        logger.log(f"[ OK ] Data file found: {file_path}\n")

        if config["type"] == "sample":
            # sample CSV file, get all columns + labels
            x, y_columns, labels = read_csv_data(file_path)
        else:
            labels = [config["label"]]
            data = np.genfromtxt(file_path)
            if config["type"] == "steady":
                x = None
                y_columns = data[1:]
            else:
                x = data[:, 0]
                y_columns = [data[:, col]]
        
    else:
        logger.log("[ ERROR ] File not found: {file_path}\n")
        
    return (x, y_columns, labels)


# RUN routine to process each test group and case
if __name__ == '__main__':

    # Construct paths
    script_dir = path.dirname(path.abspath(__file__))
    tutorials_dir = path.dirname(script_dir)
    
    # Extract current version information from file
    pmf_dir = path.dirname(tutorials_dir)
    header_file = path.join(pmf_dir, "libraries/general/PMFversion.H")
    with open(header_file, 'r') as file:
        content = file.read()
    match = re.search('const int version_ = (\S+);', content)

    if match:
        current_version = match.group(1)
    else:
        current_version = "????"

    reference_results = path.join(tutorials_dir, "referenceResults", "version.txt")
    if path.exists(reference_results):
        version_file = open(reference_results)
        reference_version = version_file.readline()
        reference_version = reference_version.replace("\n","")
        version_file.close()
    else:
        reference_version = "????"


    # Initialize logger
    logfile = path.join(tutorials_dir, 'log.plotComparison')
    logger = Logger(logfile)
    
    logger.log("========================================================")
    logger.log("                    PROCESSING DATA                     ")
    logger.log("   Reference version : " + reference_version)
    logger.log("   Current   version : " + current_version)
    logger.log("========================================================\n")

    # Create directory for figures and clean if already exists
    figures_dir = path.join(tutorials_dir, "validationFigures")

    if path.exists(figures_dir):
        shutil.rmtree(figures_dir)
    makedirs(figures_dir)
    
    for cfg in tutorials:
        solver = cfg["solver"]
        
        logger.log(f"CHECKING FOR {solver}\n")
        
        for id_case in range(len(cfg["cases"])):
            case = cfg["cases"][id_case]["case"]

            if not cfg["cases"][id_case]["result"]:
                logger.log("[ SKIP ] No validation file\n")
                continue

            if "col" in cfg["cases"][id_case]:
                col = cfg["cases"][id_case]["col"]
            else:
                col=1
            
            logger.log(f"CASE : {case}\n")
            
            # Call respective plot function based on data type    
            try:
                x, y_columns, labels = extract_XY(path.join(tutorials_dir,solver+"-tutorials"),
                                                  cfg["cases"][id_case])
                xRef, y_columnsRef, labelsRef = extract_XY(path.join(tutorials_dir,"referenceResults",solver),
                                                           cfg["cases"][id_case])

                if cfg["cases"][id_case]["type"] == "steady":
                    fig, ax = plt.subplots(figsize=(6,6), dpi=300)
                    plt.plot(y_columnsRef, y_columns, 'kx')
                    plt.plot([y_columnsRef[0],y_columnsRef[-1]],
                             [y_columnsRef[0],y_columnsRef[-1]], 'b-')
                    plt.xlabel("Reference result")
                    plt.ylabel("Current result")
                else:
                    i = 0

                    for idx, y in enumerate(y_columns):
                        fig, ax = plt.subplots(figsize=(6,6), dpi=300)
                        x_max = max(xRef)
                        y_max = max(abs(y_columnsRef[i]))
                        plt.plot(xRef/x_max, y_columnsRef[i]/y_max, label="PMF v"+reference_version, ls='-')
                        plt.plot(x/x_max, y/y_max, label="PMF v"+current_version, ls='--')
                        plt.xlabel("Dimensionless (height or time)")
                        plt.ylabel(r"Dimensionless {}".format(labels[idx]))
                        plt.legend()
                        plt.title(solver+" - "+case)
                        filename = f"{solver}_{'_'.join(case.split('/'))}_{labels[idx]}.pdf" if len(case.split('/')) > 1 else f"{solver}_{case}_{labels[idx]}.pdf"
                        plt.savefig(path.join(figures_dir,filename))
                        plt.close()
                        
                        i += 1

            except Exception as e:
                print(f"[ ERROR ] Encountered error : {e} ")
                logger.log(f"[ ERROR ] Encountered error : {e} ")

    logger.log("========================================================")
    logger.log("                        FINISHED                        ")
    logger.log("========================================================\n")
