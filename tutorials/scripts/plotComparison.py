#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 14:53:43 2024

@author: lfabregu
"""

# Imports
from __future__ import with_statement
import os, subprocess, sys, glob, time, shutil, csv
import numpy as np
import pylab as plt
plt.rc('text', usetex=False)
import re

# Import test case definitions
from tutorialsList import tutorials as testGroups

# Logger class to handle logging
class Logger:
    """Class to manage logging in a file while displaying it simultaneously."""
    def __init__(self, logfile):
        self.logfile = logfile

    def log(self, message):
        """Write message to the log file and print it."""
        with open(self.logfile, 'a') as f:
            f.write(message + "\n")

# Initialize logger
logfile = os.path.join(os.getcwd(), 'log.plot')
logger = Logger(logfile)

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
    
    pmfVersionDirectory = os.path.join('/', *res, 'solvers')
    target_file = 'headerPMF.H'
    pmfVersion = 'current_version'

    # Extract version information from file
    def extract_pmf_version(file_path):
        with open(file_path, 'r') as file:
            content = file.read()
        match = re.search(r'Current version of PMF is (\S+)', content)
        return match.group(1) if match else None

    # Verify and extract PMF version
    file_path = os.path.join(pmfVersionDirectory, target_file)
    if os.path.exists(file_path):
        pmfVersion = extract_pmf_version(file_path)
        return pmfVersion[:(len(pmfVersion)-1)]
    else:
        print('[ ERROR ] No Header Found')
        
    return pmfVersion

pmfVersion = get_pmf(os.getcwd())

# Create directory for figures and clean if already exists
figuresDir = os.path.join(os.getcwd(), "validationFigures")
if os.path.exists(figuresDir):
    shutil.rmtree(figuresDir)
os.makedirs(figuresDir)

# Function to process .gplt files and generate plots
def plot_gplt(root_dir, solver, case):
    dataDir = os.path.join(root_dir, solver, case)

    # Find .gplt file in directory
    def find_gplt_file(root_dir):
        for root, dirs, files in os.walk(root_dir):
            for file in files:
                if file.endswith('.gplt'):
                    label = file.split('.')[0].split('_')[1]
                    return (os.path.join(root, file), label)
        return None

    # Read data from .gplt file
    def read_gplt_data(file_path):
        data = []
        with open(file_path, 'r') as file:
            lines = file.readlines()
            for line in lines:
                try:
                    x, y = map(float, line.split())
                    data.append((x, y))
                except ValueError:
                    continue
        return data
    
    gpltFile = find_gplt_file(dataDir)
    if gpltFile:
        dataFile, dataLabel = gpltFile
    else:
        return("[ ERROR ] No data file found\n")

    if not dataFile:
        print("[ ERROR ] No .gplt file found in saveDatas.\n")
        logger.log("[ ERROR ] No .gplt file found in saveDatas.\n")
    else:
        print(f"[ OK ] Data file found: {dataFile}\n")
        logger.log(f"[ OK ] Data file found: {dataFile}\n")
    
        data = read_gplt_data(dataFile)
        x = np.array([point[0] for point in data])
        y = np.array([point[1] for point in data])
        
        return(x,y)
        
        plt.plot(x/max(x), y/max(abs(y)), label=r"${}$ (pmfv{})".format(dataLabel, pmfVersion))
        
        plt.xlabel("Dimensionless height")
        plt.ylabel(r"Dimensionless ${}$".format(dataLabel))
        
        plt.legend()
        if len(case.split('/')) > 1:
            strCase = '_'.join(case.split('/')).rstrip('_')
            plt.savefig(f"{solver}_{strCase}.pdf")
        else:
            plt.savefig(f"{solver}_{case}_{dataLabel}.pdf")
        
        plt.clf()

    for dir in os.listdir(os.getcwd()):
        if ".pdf" in dir:
            shutil.move(os.path.join(os.getcwd(), dir), os.path.join(figuresDir, dir))

# Function to process .csv files and generate plots
def plot_csv(root_dir, solver, case):
    dataDir = os.path.join(root_dir, solver, case)

    # Find CSV file in directory
    def find_csv_file(root_dir):
        for root, dirs, files in os.walk(root_dir):
            for file in files:
                if file.endswith('.csv'):
                    return os.path.join(root, file)
        return None

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

    csvFile = find_csv_file(dataDir)
    if csvFile:
        print(f"[ OK ] Data file found: {csvFile}\n")
        logger.log(f"[ OK ] Data file found: {csvFile}\n")

        x, y_columns, labels = read_csv_data(csvFile)
        
        
    else:
        logger.log("[ ERROR ] No .csv file found in saveDatas.\n")
        print("[ ERROR ] No .csv file found in saveDatas.\n")
        
    return (x, y_columns, labels)

# Function to process .xy files and generate plots
def plot_xy(root_dir, solver, case):
    dataDir = os.path.join(root_dir, solver, case)

    # Find .xy file in directory
    def find_xy_file(root_dir):
        for root, dirs, files in os.walk(root_dir):
            for file in files:
                if file.endswith('.xy'):
                    label = file.split('.')[0].split('_')[1]
                    return (os.path.join(root, file), label)
        return None

    # Read data from .xy file
    def read_xy_data(file_path):
        data = []
        with open(file_path, 'r') as file:
            lines = file.readlines()
            for line in lines:
                try:
                    x, y = map(float, line.split())
                    data.append((x, y))
                except ValueError:
                    continue
        return data
    
    xyFile = find_xy_file(dataDir)
    if xyFile:
        dataFile, dataLabel = xyFile
    else:
        return("[ ERROR ] No data file found\n")

    if not dataFile:
        print("[ ERROR ] No .gplt file found in saveDatas.\n")
        logger.log("[ ERROR ] No .gplt file found in saveDatas.\n")
    else:
        print(f"[ OK ] Data file found: {dataFile}\n")
        logger.log(f"[ OK ] Data file found: {dataFile}\n")
    
        data = read_xy_data(dataFile)
        x = np.array([point[0] for point in data])
        y = np.array([point[1] for point in data])
        
    return (x,y)

# RUN routine to process each test group and case
if __name__ == '__main__':
    logger.log("========================================================")
    logger.log("                    PROCESSING DATA                     ")
    logger.log("========================================================\n")

    print("========================================================")
    print("                    PROCESSING DATA                     ")
    print("========================================================\n")
    
    pmfRefVersion = ''
    
    for dir in os.listdir(os.getcwd()):
        if dir.startswith("references"):
            pmfRefVersion = dir.split('_')[1]

    for group in testGroups["tutorials"]:
        solver = group["solver"]
        k = 0
        
        print(f"CHECKING FOR {solver}\n")
        logger.log(f"CHECKING FOR {solver}\n")
        
        for case_info in group["cases"]:
            case = case_info["case"]
            dataType = group["dataTypes"][k]
            
            print(f"CASE : {case}\n")
            logger.log(f"CASE : {case}\n")
            
            # Call respective plot function based on data type    
            if dataType['dataType'] == '.csv':
                try:
                    x, y_columns, labels = plot_csv(os.path.join(os.getcwd(), 'currentResults'), solver, case)
                    xRef, y_columnsRef, labelsRef = plot_csv(os.path.join(os.getcwd(), 'referenceResults'), solver, case)
                    
                    i = 0
                    for idx, y in enumerate(y_columns):
                        plt.plot(xRef / max(xRef), y / max(abs(y_columnsRef[i])), label="reference results", ls='-')
                        plt.plot(x / max(x), y / max(abs(y)), label="current results", ls='--')
                        plt.xlabel("Dimensionless height")
                        plt.ylabel(r"Dimensionless {}".format(labels[idx]))
                        plt.legend()
            
                        filename = f"{solver}_{'_'.join(case.split('/'))}_{labels[idx]}.pdf" if len(case.split('/')) > 1 else f"{solver}_{case}_{labels[idx]}.pdf"
                        plt.savefig(filename)
                        plt.clf()
                        
                        i += 1
            
                    for pdf_file in glob.glob("*.pdf"):
                        shutil.move(pdf_file, os.path.join(figuresDir, pdf_file))
                except Exception as e:
                    print(f"[ ERROR ] Encountered error : {e} ")
                    logger.log(f"[ ERROR ] Encountered error : {e} ")
            k += 1

    logger.log("========================================================")
    logger.log("                        FINISHED                        ")
    logger.log("========================================================\n")

    print("========================================================")
    print("                        FINISHED                        ")
    print("========================================================\n")
