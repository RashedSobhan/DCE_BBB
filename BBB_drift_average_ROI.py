#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 17:19:43 2022

@author: wae16sru
"""

# code that will give average signal intensity for each ROI 
# note: when drift was calculated for all subject tthey were median drift , but once the median drift from all subjects for each ROI is achived 
# this will save the mean drift for each ROI (in two columns: 1. overall drift includoing the negative value, and 2. positive drift including only the poitive drift )
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from numpy import trapz
from sklearn.linear_model import LinearRegression

model = LinearRegression()

df_En = pd.read_csv('/gpfs/home/wae16sru/BBB_sample_data/BBB_Michael_Code_modification/BBB_manuscript_results/BBB_drift_results/BBB_drift_sigEn_allsub_on_T1_std.csv')

df_In = pd.read_csv('/gpfs/home/wae16sru/BBB_sample_data/BBB_Michael_Code_modification/BBB_manuscript_results/BBB_drift_results/BBB_drift_sigIn_allsub_on_T1_std.csv')

df_In_norm= pd.read_csv('/gpfs/home/wae16sru/BBB_sample_data/BBB_Michael_Code_modification/BBB_manuscript_results/BBB_drift_results/BBB_drift_sigIn_norm_allsub_on_T1_std.csv')


df_En.columns= ["ROIS", "t0", "t1", "t2", "t3", "t4", "t5","t6","t7", \
                                   "t8","t9", "t10", "t11", "t12", "t13", "t14", "t15",\
                                       "t16", "t17", "t18", "t19"]
df_In.columns= ["ROIS", "t0", "t1", "t2", "t3", "t4", "t5","t6","t7", \
                                   "t8","t9", "t10", "t11", "t12", "t13", "t14", "t15",\
                                       "t16", "t17", "t18", "t19"]

df_In_norm.columns= ["ROIS", "t0", "t1", "t2", "t3", "t4", "t5","t6","t7", \
                                   "t8","t9", "t10", "t11", "t12", "t13", "t14", "t15",\
                                       "t16", "t17", "t18", "t19"]


results_dir="/gpfs/home/wae16sru/BBB_sample_data/BBB_Michael_Code_modification/BBB_manuscript_results/BBB_drift_results"
# TIME SERIES SIGNAL ENHANCEMENT CALCULATIONS
#==============================================
# Calculate timings of each acquisition time point. ** USER INPUT **
acq_len = 591.5  # Acquisition duration (s) (subtracted 31.5 so pt 0 = time 0)
nPts = len(df_En.columns)-1
tVec = np.arange( 0, nPts ) * acq_len / nPts

print(df_En.head())
print(df_In.head())

# sigEn_tissue_drift=df.groupby(['ROIS']).mean()
sigEn_tissue_drift= df_En.groupby(['ROIS']).median(0)
sigIn_tissue_drift= df_In.groupby(['ROIS']).median(0)
sigIn_norm_tissue_drift= df_In_norm.groupby(['ROIS']).median(0)


sigEn_tissue_drift_mean= df_En.groupby(['ROIS']).mean(0)
sigIn_tissue_drift_mean= df_In.groupby(['ROIS']).mean(0)
sigIn_norm_tissue_drift_mean= df_In_norm.groupby(['ROIS']).mean(0)


print (sigEn_tissue_drift)

print(sigEn_tissue_drift.loc[['WM'],:])


print(sigEn_tissue_drift[0:1])


print(sigEn_tissue_drift.loc['WM', sigEn_tissue_drift.columns[0:20]])

sigEn_WM_drift_mean= sigEn_tissue_drift_mean.loc['WM', sigEn_tissue_drift_mean.columns[0:20]].to_numpy(dtype=float)
sigEn_GM_drift_mean= sigEn_tissue_drift_mean.loc['GM', sigEn_tissue_drift_mean.columns[0:20]].to_numpy(dtype=float)
sigEn_CSF_drift_mean= sigEn_tissue_drift_mean.loc['CSF', sigEn_tissue_drift_mean.columns[0:20]].to_numpy(dtype=float)
sigEn_Hippo_drift_mean= sigEn_tissue_drift_mean.loc['Hippo', sigEn_tissue_drift_mean.columns[0:20]].to_numpy(dtype=float)
sigEn_Thalamus_drift_mean= sigEn_tissue_drift_mean.loc['Thalamus', sigEn_tissue_drift_mean.columns[0:20]].to_numpy(dtype=float)
sigEn_Amygdala_drift_mean= sigEn_tissue_drift_mean.loc['Amygdala', sigEn_tissue_drift_mean.columns[0:20]].to_numpy(dtype=float)

sigIn_WM_drift_mean= sigIn_tissue_drift_mean.loc['WM', sigEn_tissue_drift_mean.columns[0:20]].to_numpy(dtype=float)
sigIn_GM_drift_mean = sigIn_tissue_drift_mean.loc['GM', sigEn_tissue_drift_mean.columns[0:20]].to_numpy(dtype=float)
sigIn_CSF_drift_mean = sigIn_tissue_drift_mean.loc['CSF', sigEn_tissue_drift_mean.columns[0:20]].to_numpy(dtype=float)
sigIn_Hippo_drift_mean = sigIn_tissue_drift_mean.loc['Hippo', sigEn_tissue_drift_mean.columns[0:20]].to_numpy(dtype=float)
sigIn_Thalamus_drift_mean = sigIn_tissue_drift_mean.loc['Thalamus', sigEn_tissue_drift_mean.columns[0:20]].to_numpy(dtype=float)
sigIn_Amygdala_drift_mean = sigIn_tissue_drift_mean.loc['Amygdala', sigEn_tissue_drift_mean.columns[0:20]].to_numpy(dtype=float)

sigIn_WM_norm_drift_mean = sigIn_norm_tissue_drift_mean.loc['WM', sigEn_tissue_drift_mean.columns[0:20]].to_numpy(dtype=float)
sigin_GM_norm_drift_mean = sigIn_norm_tissue_drift_mean.loc['GM', sigEn_tissue_drift_mean.columns[0:20]].to_numpy(dtype=float)
sigIn_CSF_norm_drift_mean = sigIn_norm_tissue_drift_mean.loc['CSF', sigEn_tissue_drift_mean.columns[0:20]].to_numpy(dtype=float)
sigIn_Hippo_norm_drift_mean = sigIn_norm_tissue_drift_mean.loc['Hippo', sigEn_tissue_drift_mean.columns[0:20]].to_numpy(dtype=float)
sigIn_Thalamus_norm_drift_mean = sigIn_norm_tissue_drift_mean.loc['Thalamus', sigEn_tissue_drift_mean.columns[0:20]].to_numpy(dtype=float)
sigIn_Amygdala_norm_drift_mean = sigIn_norm_tissue_drift_mean.loc['Amygdala', sigEn_tissue_drift_mean.columns[0:20]].to_numpy(dtype=float)


#same thing to do for intensity 

#print (sigIn_tissue_drift)

#rint(sigIn_tissue_drift.loc[['WM'],:])


#print(sigIn_tissue_drift[0:1])


#print(sigIn_tissue_drift.loc['WM', sigIn_tissue_drift.columns[0:20]])

sigEn_WM_drift= sigEn_tissue_drift.loc['WM', sigIn_tissue_drift.columns[0:20]].to_numpy(dtype=float)
sigEn_GM_drift= sigEn_tissue_drift.loc['GM', sigIn_tissue_drift.columns[0:20]].to_numpy(dtype=float)
sigEn_CSF_drift= sigEn_tissue_drift.loc['CSF', sigIn_tissue_drift.columns[0:20]].to_numpy(dtype=float)
sigEn_Hippo_drift= sigEn_tissue_drift.loc['Hippo', sigIn_tissue_drift.columns[0:20]].to_numpy(dtype=float)
sigEn_Thalamus_drift= sigEn_tissue_drift.loc['Thalamus', sigIn_tissue_drift.columns[0:20]].to_numpy(dtype=float)
sigEn_Amygdala_drift= sigEn_tissue_drift.loc['Amygdala', sigIn_tissue_drift.columns[0:20]].to_numpy(dtype=float)

sigIn_WM_drift= sigIn_tissue_drift.loc['WM', sigIn_tissue_drift.columns[0:20]].to_numpy(dtype=float)
sigIn_GM_drift= sigIn_tissue_drift.loc['GM', sigIn_tissue_drift.columns[0:20]].to_numpy(dtype=float)
sigIn_CSF_drift = sigIn_tissue_drift.loc['CSF', sigIn_tissue_drift.columns[0:20]].to_numpy(dtype=float)
sigIn_Hippo_drift = sigIn_tissue_drift.loc['Hippo', sigIn_tissue_drift.columns[0:20]].to_numpy(dtype=float)
sigIn_Thalamus_drift = sigIn_tissue_drift.loc['Thalamus', sigIn_tissue_drift.columns[0:20]].to_numpy(dtype=float)
sigIn_Amygdala_drift = sigIn_tissue_drift.loc['Amygdala', sigIn_tissue_drift.columns[0:20]].to_numpy(dtype=float)

sigIn_norm_WM_drift= sigIn_norm_tissue_drift.loc['WM', sigIn_norm_tissue_drift.columns[0:20]].to_numpy(dtype=float)
sigIn_norm_GM_drift= sigIn_norm_tissue_drift.loc['GM', sigIn_norm_tissue_drift.columns[0:20]].to_numpy(dtype=float)
sigIn_norm_CSF_drift = sigIn_norm_tissue_drift.loc['CSF', sigIn_norm_tissue_drift.columns[0:20]].to_numpy(dtype=float)
sigIn_norm_Hippo_drift = sigIn_norm_tissue_drift.loc['Hippo', sigIn_norm_tissue_drift.columns[0:20]].to_numpy(dtype=float)
sigIn_norm_Thalamus_drift = sigIn_norm_tissue_drift.loc['Thalamus', sigIn_norm_tissue_drift.columns[0:20]].to_numpy(dtype=float)
sigIn_norm_Amygdala_drift = sigIn_norm_tissue_drift.loc['Amygdala', sigIn_norm_tissue_drift.columns[0:20]].to_numpy(dtype=float)


# Plotting code - shows signal intensity time-curves for WM, GM, CSF, sigEn_SS SS
fig = plt.figure( )
fig.suptitle( 'Median Signal drift per ROI ', fontsize = 14 )
ax = plt.gca( )
plt1, = plt.plot( tVec, sigEn_WM_drift * 100, '#56B4E9', label = 'White Matter', \
linewidth = 2 )
plt2, = plt.plot( tVec, sigEn_GM_drift * 100, '#0072B2', label = 'Grey Matter' , \
linewidth = 2 )
plt3, = plt.plot( tVec, sigEn_Hippo_drift * 100, '#CC79A7', label = 'Hippocampus' , \
linewidth = 2 )
plt4, = plt.plot( tVec, sigEn_Thalamus_drift * 100, '#00FFFF', label = 'Thalamus' , \
linewidth = 2 )
plt5, = plt.plot( tVec, sigEn_Amygdala_drift * 100, '#8A2BE2', label = 'Amygdala' , \
linewidth = 2 )
plt6, = plt.plot( tVec, sigEn_CSF_drift * 100, '#8B7355', label ='CSF' , \
linewidth = 2 )
#????more regions to add here 

ax.set_xlabel( 'Time (s)' )
ax.set_ylabel( ' Signal Enhancement (drift) (%)' )
ax.grid(None)
# Put SS on similar scale

plt.legend( handles = [ plt1, plt2, plt3, plt4, plt5, plt6], loc = 'lower center', mode = "expand", ncol = 3,frameon = False ) 
#????more regions to add here 
plt.show( )
plt.savefig(results_dir+'/BBB_ROI_medianEn_drift_plot_T1_std.png', format = 'png', bbox_inches = 'tight', dpi = 500)



#not doing this stage for intensity (maybe later)
#doing this for normalised intensity 

# Plotting code - shows signal intensity time-curves for WM, GM, CSF, sigEn_SS SS
fig = plt.figure( )
fig.suptitle( 'Median Signal drift per ROI', fontsize = 14 )
ax = plt.gca( )
plt1, = plt.plot( tVec, sigIn_norm_WM_drift * 100, '#56B4E9', label = 'White Matter', \
linewidth = 2 )
plt2, = plt.plot( tVec, sigIn_norm_GM_drift* 100, '#0072B2', label = 'Grey Matter' , \
linewidth = 2 )
plt3, = plt.plot( tVec, sigIn_norm_Hippo_drift* 100, '#CC79A7', label = 'Hippocampus' , \
linewidth = 2 )
plt4, = plt.plot( tVec, sigIn_norm_Thalamus_drift* 100, '#00FFFF', label = 'Thalamus' , \
linewidth = 2 )
plt5, = plt.plot( tVec, sigIn_norm_Amygdala_drift* 100, '#8A2BE2', label = 'Amygdala' , \
linewidth = 2 )
plt6, = plt.plot( tVec, sigIn_norm_CSF_drift* 100, '#8B7355', label ='CSF' , \
linewidth = 2 )
#????more regions to add here 

ax.set_xlabel( 'Time (s)' )
ax.set_ylabel( ' Normalised Signal intensity (%)' )
ax.grid(None)


plt.legend( handles = [ plt1, plt2, plt3, plt4, plt5, plt6], loc = 'lower center', mode = "expand", ncol = 3,frameon = False ) 
#????more regions to add here 
plt.show( )
plt.savefig(results_dir+'/BBB_ROI_mediansigIn_norm_drift_plot_T1_std.png', format = 'png', bbox_inches = 'tight', dpi = 500)

drift_En_allsub= pd.read_csv(results_dir+'/BBB_drift_polyfit_sigEn_on_T1_std.csv')
drift_En_allsub.columns=["ROIS", "Drift_Enchancement_POLYFIT"]
drift_ROI_average= drift_En_allsub.groupby(['ROIS']).agg([np.average])

drift_In_norm_allsub= pd.read_csv(results_dir+'/BBB_drift_polyfit_sigIn_norm_on_T1_std.csv')
drift_In_norm_allsub.columns=["ROIS", "Drift_normalised_signal_POLYFIT"]
drift_In_ROI_average= drift_In_norm_allsub.groupby(['ROIS']).agg([np.average])

#driftI_average.to_csv('/gpfs/home/wae16sru/BBB_sample_data/BBB_Michael_Code_modification/drift_ROI_mean_polyfit.csv', mode='a', header= True)

one_signal_unit_mean=[1/sigIn_Amygdala_drift_mean[0],1/sigIn_CSF_drift_mean[0],
1/sigIn_GM_drift_mean[0],1/sigIn_Hippo_drift_mean[0],1/sigIn_Thalamus_drift_mean[0],1/sigIn_WM_drift_mean[0]]

one_signal_unit_median=[1/sigIn_Amygdala_drift[0],1/sigIn_CSF_drift[0],
1/sigIn_GM_drift[0],1/sigIn_Hippo_drift[0],1/sigIn_Thalamus_drift[0],1/sigIn_WM_drift[0]]

drift_ROI_average['Drift_normalised_signal_POLYFIT']=drift_In_ROI_average.loc[:,"Drift_normalised_signal_POLYFIT"]
#drift_ROI_averge['Drift_normalised signal_polyfit']=drift_In_norm_allsub[1:][1]
drift_ROI_average['One_Signal_Unit_Mean'] = one_signal_unit_mean
drift_ROI_average['One_Signal_Unit_Median'] = one_signal_unit_median

print (drift_ROI_average)
drift_ROI_average.to_csv(results_dir+'/BBB_ROI_drift_mean_polyfit.csv', mode='a', header= True)