#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 11:31:44 2022

@author: phorgue
"""

import numpy as np
import matplotlib.pyplot as plt

data1 = np.genfromtxt('postProcessing/sampleDict/864/acrossFlow_C_CMatrix_theta_thetaMatrix.csv', delimiter=',')
data1 = data1[1:,:]
data4 = np.genfromtxt('postProcessing/sampleDict/3456/acrossFlow_C_CMatrix_theta_thetaMatrix.csv', delimiter=',')
data4 = data4[1:,:]
data8 = np.genfromtxt('postProcessing/sampleDict/6912/acrossFlow_C_CMatrix_theta_thetaMatrix.csv', delimiter=',')
data8 = data8[1:,:]

data_metis_flow_1f = np.genfromtxt('metis_results/fracture_saturation_t864.txt')
data_metis_flow_4f = np.genfromtxt('metis_results/fracture_saturation_t3456.txt')
data_metis_flow_8f = np.genfromtxt('metis_results/fracture_saturation_t6912.txt')
data_metis_flow_1m = np.genfromtxt('metis_results/matrix_saturation_t864.txt')
data_metis_flow_4m = np.genfromtxt('metis_results/matrix_saturation_t3456.txt')
data_metis_flow_8m = np.genfromtxt('metis_results/matrix_saturation_t6912.txt')

data_metis_tracer_1f = np.genfromtxt('metis_results/fracture_tracer_t864.txt')
data_metis_tracer_4f = np.genfromtxt('metis_results/fracture_tracer_t3456.txt')
data_metis_tracer_8f = np.genfromtxt('metis_results/fracture_tracer_t6912.txt')
data_metis_tracer_1m = np.genfromtxt('metis_results/matrix_tracer_t864.txt')
data_metis_tracer_4m = np.genfromtxt('metis_results/matrix_tracer_t3456.txt')
data_metis_tracer_8m = np.genfromtxt('metis_results/matrix_tracer_t6912.txt')

fig1, ax1 = plt.subplots(figsize=(8,8))
ax1.set_xlabel('Water content [%]')
ax1.set_ylabel('Height [m]')
plt.xlim([0,0.025])
plt.ylim([0,0.45])
plt.plot(data1[:,3], data1[:,0], 'b-', label='PMF T=0.01 day')
plt.plot(data_metis_flow_1f[:,1]*0.025, data_metis_flow_1f[:,0], 'b--', label='METIS T=0.01 day')
plt.plot(data4[:,3], data4[:,0], 'r-', label='PMF T=0.04 day')
plt.plot(data_metis_flow_4f[:,1]*0.025, data_metis_flow_4f[:,0], 'r--', label='METIS T=0.04 day')
plt.plot(data8[:,3], data8[:,0], 'g-', label='PMF T=0.08 day')
plt.plot(data_metis_flow_8f[:,1]*0.025, data_metis_flow_8f[:,0], 'g--', label='METIS T=0.08 day')
ax1.legend()
fig1.savefig("results_flow_fracture.pdf")

fig2, ax2 = plt.subplots(figsize=(8,8))
ax2.set_xlabel('Water content [%]')
ax2.set_ylabel('Height [m]')
plt.xlim([0.25,0.40])
plt.ylim([0,0.45])
plt.plot(data1[:,4], data1[:,0], 'b-', label='PMF T=0.01 jour')
plt.plot(data_metis_flow_1m[:,1]*0.475, data_metis_flow_1m[:,0], 'b--', label='METIS T=0.01 day')
plt.plot(data4[:,4], data4[:,0], 'r-', label='PMF T=0.04 jour')
plt.plot(data_metis_flow_4m[:,1]*0.475, data_metis_flow_4m[:,0], 'r--', label='METIS T=0.04 day')
plt.plot(data8[:,4], data8[:,0], 'g-', label='PMF T=0.08 jour')
plt.plot(data_metis_flow_8m[:,1]*0.475, data_metis_flow_8m[:,0], 'g--', label='METIS T=0.08 day')
ax2.legend()
fig2.savefig("results_flow_matrix.pdf")

fig3, ax3 = plt.subplots(figsize=(8,8))
ax3.set_xlabel('Concentration [g/L]')
ax3.set_ylabel('Height [m]')
plt.xlim([0,1.2])
plt.ylim([0,0.45])
plt.plot(data1[:,1], data1[:,0], 'b-', label='PMF T=0.01 day')
plt.plot(data_metis_tracer_1f[:,1], data_metis_tracer_1f[:,0], 'b--', label='METIS T=0.01 day')
plt.plot(data4[:,1], data4[:,0], 'r-', label='PMF T=0.04 day')
plt.plot(data_metis_tracer_4f[:,1], data_metis_tracer_4f[:,0], 'r--', label='METIS T=0.04 day')
plt.plot(data8[:,1], data8[:,0], 'g-', label='PMF T=0.08 day')
plt.plot(data_metis_tracer_8f[:,1], data_metis_tracer_8f[:,0], 'g--', label='METIS T=0.08 day')
ax3.legend()
fig3.savefig("results_tracer_fracture.pdf")

fig4, ax4 = plt.subplots(figsize=(8,8))
ax4.set_xlabel('Concentration [g/L]')
ax4.set_ylabel('Height [m]')
plt.xlim([0,1.2])
plt.ylim([0,0.45])
plt.plot(data1[:,2], data1[:,0], 'b-', label='PMF T=0.01 jour')
plt.plot(data_metis_tracer_1m[:,1], data_metis_tracer_1m[:,0], 'b--', label='METIS T=0.01 day')
plt.plot(data4[:,2], data4[:,0], 'r-', label='PMF T=0.04 jour')
plt.plot(data_metis_tracer_4m[:,1], data_metis_tracer_4m[:,0], 'r--', label='METIS T=0.04 day')
plt.plot(data8[:,2], data8[:,0], 'g-', label='PMF T=0.08 jour')
plt.plot(data_metis_tracer_8m[:,1], data_metis_tracer_8m[:,0], 'g--', label='METIS T=0.08 day')
ax4.legend()
fig4.savefig("results_tracer_matrix.pdf")
