#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 11:54:38 2023

@author: wae16sru
"""

#this file will crate colopmap for t1_HIFI, kappa_HIFI and s0_HIFI

import os

sub_ID= input("enter subject ID: ")
sub_ID_split=sub_ID.split("_")
sub_ID_short=int(sub_ID_split[0])

dir= "/gpfs/home/wae16sru/BBB_sample_data/" + sub_ID + "/T1-mapping_results_RS"
os.chdir(str(dir))
results_dir= "/gpfs/home/wae16sru/BBB_sample_data/BBB_Michael_Code_modification/BBB_manuscript_results/BBB_DCE_results"



import sys
import numpy as np
import matplotlib.pyplot as plt


import math
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt


plt.style.use( 'ggplot' )
plt.ion()

t1_HIFI = nib.load(dir+'/t1_HIFI.nii')
t1_HIFI_data = t1_HIFI.get_fdata( )
hdr_T1img = t1_HIFI.header


kappa_HIFI = nib.load(dir+'/k_fa_HIFI.nii')
kappa_HIFI_data = kappa_HIFI.get_fdata( )
hdr_kappaimg = kappa_HIFI.header

s0_HIFI = nib.load(dir+'/s0_HIFI.nii')
s0_HIFI_data = s0_HIFI.get_fdata( )
hdr_s0img = s0_HIFI.header



img_T1=t1_HIFI_data[:,:,103]

img_kappa=kappa_HIFI_data[:,:,103]
img_s0=s0_HIFI_data[:,:,103]


#making colorbar directions: https://matplotlib.org/stable/gallery/color/colorbar_basics.html

fig, (ax1, ax2, ax3) = plt.subplots(figsize=(13,3), ncols=3)


#plot the t1_data
img_T1_colormap = ax1.imshow(img_T1, cmap='magma', vmin=0, vmax=10, interpolation='none', aspect='auto')
fig.colorbar(img_T1_colormap, ax=ax1)

#plot the kappa_map

img_kappa_colormap = ax2.imshow(img_kappa, cmap='plasma', vmin=0, vmax=3, interpolation='none',aspect='auto')
fig.colorbar(img_kappa_colormap, ax=ax2)


img_s0_colormap = ax3.imshow(img_s0, cmap='Greys', vmin=800, vmax=100000, interpolation='none',aspect='auto')
fig.colorbar(img_s0_colormap, ax=ax3)