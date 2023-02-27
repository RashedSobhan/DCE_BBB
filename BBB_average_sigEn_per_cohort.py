#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 17:12:50 2023

@author: wae16sru
"""

# this code is to create average signal enhancement from each cohort 

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


results_dir= "/gpfs/home/wae16sru/BBB_sample_data/BBB_Michael_Code_modification/BBB_manuscript_results/BBB_DCE_results"
df_sigEn_SMI_control = pd.read_csv('/gpfs/home/wae16sru/BBB_sample_data/BBB_Michael_Code_modification/BBB_manuscript_results/BBB_DCE_results/all_sigEn_SMI_control.csv')
df_sigEn_SMI_active =  pd.read_csv('/gpfs/home/wae16sru/BBB_sample_data/BBB_Michael_Code_modification/BBB_manuscript_results/BBB_DCE_results/all_sigEn_SMI_active.csv')
df_sigEn_MCI_control = pd.read_csv('/gpfs/home/wae16sru/BBB_sample_data/BBB_Michael_Code_modification/BBB_manuscript_results/BBB_DCE_results/all_sigEn_MCI_control.csv')



df_sigEn_SMI_control.columns= ["ROIS", "t0", "t1", "t2", "t3", "t4", "t5","t6","t7", \
                                   "t8","t9", "t10", "t11", "t12", "t13", "t14", "t15",\
                                       "t16", "t17", "t18", "t19"]

df_sigEn_SMI_active.columns= ["ROIS", "t0", "t1", "t2", "t3", "t4", "t5","t6","t7", \
                                   "t8","t9", "t10", "t11", "t12", "t13", "t14", "t15",\
                                       "t16", "t17", "t18", "t19"]
    
df_sigEn_MCI_control.columns= ["ROIS", "t0", "t1", "t2", "t3", "t4", "t5","t6","t7", \
                                   "t8","t9", "t10", "t11", "t12", "t13", "t14", "t15",\
                                       "t16", "t17", "t18", "t19"]
    
median_sigEn_SMI_control= df_sigEn_SMI_control.groupby(['ROIS']).median(0)
median_sigEn_SMI_active= df_sigEn_SMI_active.groupby(['ROIS']).median(0)
median_sigEn_MCI_control= df_sigEn_MCI_control.groupby(['ROIS']).median(0)

# TIME SERIES SIGNAL ENHANCEMENT CALCULATIONS
#==============================================
# Calculate timings of each acquisition time point. ** USER INPUT **
acq_len = 591.5  # Acquisition duration (s) (subtracted 31.5 so pt 0 = time 0)
nPts = len(df_sigEn_MCI_control.columns)-1
tVec = np.arange( 0, nPts ) * acq_len / nPts

sigEn_WM_SMI_control= median_sigEn_SMI_control.loc['WM', median_sigEn_SMI_control.columns[0:20]].to_numpy(dtype=float)
sigEn_GM_SMI_control= median_sigEn_SMI_control.loc['GM', median_sigEn_SMI_control.columns[0:20]].to_numpy(dtype=float)
sigEn_CSF_SMI_control= median_sigEn_SMI_control.loc['CSF', median_sigEn_SMI_control.columns[0:20]].to_numpy(dtype=float)
sigEn_Hippo_SMI_control= median_sigEn_SMI_control.loc['Hippo', median_sigEn_SMI_control.columns[0:20]].to_numpy(dtype=float)
sigEn_Thalamus_SMI_control= median_sigEn_SMI_control.loc['Thalamus', median_sigEn_SMI_control.columns[0:20]].to_numpy(dtype=float)
sigEn_Amygdala_SMI_control= median_sigEn_SMI_control.loc['Amygdala', median_sigEn_SMI_control.columns[0:20]].to_numpy(dtype=float)


sigEn_WM_SMI_active= median_sigEn_SMI_active.loc['WM', median_sigEn_SMI_active.columns[0:20]].to_numpy(dtype=float)
sigEn_GM_SMI_active= median_sigEn_SMI_active.loc['GM', median_sigEn_SMI_active.columns[0:20]].to_numpy(dtype=float)
sigEn_CSF_SMI_active= median_sigEn_SMI_active.loc['CSF', median_sigEn_SMI_active.columns[0:20]].to_numpy(dtype=float)
sigEn_Hippo_SMI_active= median_sigEn_SMI_active.loc['Hippo', median_sigEn_SMI_active.columns[0:20]].to_numpy(dtype=float)
sigEn_Thalamus_SMI_active= median_sigEn_SMI_active.loc['Thalamus', median_sigEn_SMI_active.columns[0:20]].to_numpy(dtype=float)
sigEn_Amygdala_SMI_active= median_sigEn_SMI_active.loc['Amygdala', median_sigEn_SMI_active.columns[0:20]].to_numpy(dtype=float)

sigEn_WM_MCI_control= median_sigEn_MCI_control.loc['WM', median_sigEn_MCI_control.columns[0:20]].to_numpy(dtype=float)
sigEn_GM_MCI_control= median_sigEn_MCI_control.loc['GM', median_sigEn_MCI_control.columns[0:20]].to_numpy(dtype=float)
sigEn_CSF_MCI_control= median_sigEn_MCI_control.loc['CSF', median_sigEn_MCI_control.columns[0:20]].to_numpy(dtype=float)
sigEn_Hippo_MCI_control= median_sigEn_MCI_control.loc['Hippo', median_sigEn_MCI_control.columns[0:20]].to_numpy(dtype=float)
sigEn_Thalamus_MCI_control= median_sigEn_MCI_control.loc['Thalamus', median_sigEn_MCI_control.columns[0:20]].to_numpy(dtype=float)
sigEn_Amygdala_MCI_control= median_sigEn_MCI_control.loc['Amygdala', median_sigEn_MCI_control.columns[0:20]].to_numpy(dtype=float)


# Plotting code - shows signal intensity time-curves for WM, GM, CSF, sigEn_SS SS
fig = plt.figure( )
fig.suptitle( 'Median Signal drift per ROI (SMI-control)', fontsize = 14 )
ax = plt.gca( )
plt1, = plt.plot( tVec, sigEn_WM_SMI_control, '#56B4E9', label = 'White Matter', \
linewidth = 2 )
plt2, = plt.plot( tVec, sigEn_GM_SMI_control, '#0072B2', label = 'Grey Matter' , \
linewidth = 2 )
plt3, = plt.plot( tVec, sigEn_Hippo_SMI_control, '#CC79A7', label = 'Hippocampus' , \
linewidth = 2 )
plt4, = plt.plot( tVec, sigEn_Thalamus_SMI_control, '#00FFFF', label = 'Thalamus' , \
linewidth = 2 )
plt5, = plt.plot( tVec, sigEn_Amygdala_SMI_control, '#8A2BE2', label = 'Amygdala' , \
linewidth = 2 )
plt6, = plt.plot( tVec, sigEn_CSF_SMI_control, '#8B7355', label ='CSF' , \
linewidth = 2 )
#????more regions to add here 

ax.set_xlabel( 'Time (s)' )
ax.set_ylabel( ' Signal Enhancement (%)' )
ax.grid(None)
# Put SS on similar scale

plt.legend( handles = [ plt1, plt2, plt3, plt4, plt5, plt6], loc = 'lower center', mode = "expand", ncol = 3,frameon = False ) 
#????more regions to add here 
plt.show( )
plt.savefig(results_dir+'/BBB_median_sigEn_SMI_control.png', format = 'png', bbox_inches = 'tight', dpi = 500)


fig = plt.figure( )
fig.suptitle( 'Median Signal drift per ROI (SMI-active)', fontsize = 14 )
ax = plt.gca( )
plt1, = plt.plot( tVec, sigEn_WM_SMI_active, '#56B4E9', label = 'White Matter', \
linewidth = 2 )
plt2, = plt.plot( tVec, sigEn_GM_SMI_active, '#0072B2', label = 'Grey Matter' , \
linewidth = 2 )
plt3, = plt.plot( tVec, sigEn_Hippo_SMI_active, '#CC79A7', label = 'Hippocampus' , \
linewidth = 2 )
plt4, = plt.plot( tVec, sigEn_Thalamus_SMI_active, '#00FFFF', label = 'Thalamus' , \
linewidth = 2 )
plt5, = plt.plot( tVec, sigEn_Amygdala_SMI_active, '#8A2BE2', label = 'Amygdala' , \
linewidth = 2 )
plt6, = plt.plot( tVec, sigEn_CSF_SMI_active, '#8B7355', label ='CSF' , \
linewidth = 2 )
#????more regions to add here 

ax.set_xlabel( 'Time (s)' )
ax.set_ylabel( ' Signal Enhancement (%)' )
ax.grid(None)
# Put SS on similar scale

plt.legend( handles = [ plt1, plt2, plt3, plt4, plt5, plt6], loc = 'lower center', mode = "expand", ncol = 3,frameon = False ) 
#????more regions to add here 
plt.show( )
plt.savefig(results_dir+'/BBB_median_sigEn_SMI_active.png', format = 'png', bbox_inches = 'tight', dpi = 500)

# Plotting code - shows signal intensity time-curves for WM, GM, CSF, sigEn_SS SS
fig = plt.figure( )
fig.suptitle( 'Median Signal drift per ROI (MCI-control)', fontsize = 14 )
ax = plt.gca( )
plt1, = plt.plot( tVec, sigEn_WM_MCI_control, '#56B4E9', label = 'White Matter', \
linewidth = 2 )
plt2, = plt.plot( tVec, sigEn_GM_MCI_control, '#0072B2', label = 'Grey Matter' , \
linewidth = 2 )
plt3, = plt.plot( tVec, sigEn_Hippo_MCI_control, '#CC79A7', label = 'Hippocampus' , \
linewidth = 2 )
plt4, = plt.plot( tVec, sigEn_Thalamus_MCI_control, '#00FFFF', label = 'Thalamus' , \
linewidth = 2 )
plt5, = plt.plot( tVec, sigEn_Amygdala_MCI_control, '#8A2BE2', label = 'Amygdala' , \
linewidth = 2 )
plt6, = plt.plot( tVec, sigEn_CSF_MCI_control, '#8B7355', label ='CSF' , \
linewidth = 2 )
#????more regions to add here 

ax.set_xlabel( 'Time (s)' )
ax.set_ylabel( ' Signal Enhancement (%)' )
ax.grid(None)
# Put SS on similar scale

plt.legend( handles = [ plt1, plt2, plt3, plt4, plt5, plt6], loc = 'lower center', mode = "expand", ncol = 3,frameon = False ) 
#????more regions to add here 
plt.show( )
plt.savefig(results_dir+'/BBB_median_sigEn_MCI_control.png', format = 'png', bbox_inches = 'tight', dpi = 500)

