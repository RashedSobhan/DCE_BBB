#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 11:44:54 2022

@author: wae16sru
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 22:46:51 2022

@author: wae16sru (Rashed Sobhan) 

"""
#inspired from Michale's code from SEPAL : https://github.com/mjt320/SEPAL


# this code does t1-mapping according to SEPAL 
 
import os

sub_ID= input("enter subject ID: ")
sub_ID_split=sub_ID.split("_")
sub_ID_short=int(sub_ID_split[0])

dir= "/gpfs/home/wae16sru/BBB_sample_data/" + sub_ID
os.chdir(str(dir))


import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('/gpfs/home/wae16sru/BBB_sample_data/BBB_Michael_Code_modification/SEPAL-master/src')
import t1_fit,dce_fit, relaxivity, signal_models, water_ex_models, aifs, pk_models

FA5_dir=dir+"/products_RS/DESPOT1_5_reg.nii.gz"
FA10_dir=dir+"/products_RS/DESPOT1_10_reg.nii.gz"
FA15_dir=dir+"/products_RS/DESPOT1_15_reg.nii.gz"
TI_dir=dir+"/products_RS/DESPOT1_IR5_reg.nii.gz"


esp = np.array([5.564e-3, 5.568e-3, 5.568e-3, 5.568e-3]) # echo spacing (IR-SPGR) or TR (SPGR) ?????
# esp = np.array([1.5e-3, 5.568e-3, 5.568e-3, 5.568e-3]) 
ti = np.array([0.45, np.nan, np.nan, np.nan]) # delay after inversion pulse
n = np.array([96, np.nan, np.nan, np.nan]) # number of readout pulses (IR-SPGR only)
b = np.array([5, 5, 10, 15]) # excitation flip angle
td = np.array([0, np.nan, np.nan, np.nan]) # delay between end of readout train and next inversion pulse (IR-SPGR only)
centre = np.array([0, np.nan, np.nan, np.nan]) # time when centre of k-space is acquired (expressed as fraction of readout pulse train length; IR-SPGR only)

#images = [os.path.join('.', 'T1_data', img) for img in ['TI_168ms.nii.gz', 'TI_1068ms.nii.gz', 'FA2.nii.gz', 'FA5.nii.gz', 'FA12.nii.gz']]

hifi_fitter = t1_fit.HIFI(esp, ti, n, b, td, centre)

images = [os.path.join('.', 'dir', img) for img in [TI_dir, FA5_dir, FA10_dir, FA15_dir]]


s0, t1, k_fa, s_fit = hifi_fitter.proc_image(images, 
                                             mask=dir+'/brain_mask_FINAL.nii.gz',
                                             threshold=50,
                                             dir=dir+'/T1-mapping_results_RS',
                                             suffix='_HIFI',
                                             n_procs=4);
