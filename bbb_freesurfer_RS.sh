#!/bin/bash
bbb_pipe(){
#=========================================================================
# TITLE: "BBB Segementation Pipeline, Version 0.1" (bbb_fastsurfer.sh)
# PURPOSE: Script for obtaining ROI from which signal enhancement will be extracted. 
# this code Performs segmentation of T1-weighted images in FreeSurfer. 
#
# INPUT: 
# AUTHORS: Donnie Cameron & Rashed Sobhan
# CREATED: 2022/05/05
# LAST UPDATED:  01/12/2022
#
# REQUIRED PACKAGES: FreeSurfer 6.0.0, FSL 6.0.0
#=========================================================================

# Before commencing, convert essential DICOM files to Nifti format using 
# dcm2niix (faster than dcm2nii and it includes GE sequence names in 
# filenames). dcm2niix can be downloaded as part of MRIcron

# DCam - Check dcm2niix-converted TRACC data and adapt as needed...
#        Also need to allow option for manual CBF factor calc   

# Error if function called with no arguments...
local val=${1:?Must provide at least one argument}
echo "$val"



echo "bbb_pipe(): Dataset (\$@) passed - \"$@\""

subj_id="$@"
module add freesurfer
export SUBJECTS_DIR=~/BBB_sample_data/$subj_id
. /gpfs/software/ada/freesurfer/6.0/FreeSurferEnv.sh
module add fsl



# Create directories for FSL products and QC images
if [ ! -d $SUBJECTS_DIR/products2 ]
then
    mkdir $SUBJECTS_DIR/products2
fi
if [ ! -d $SUBJECTS_DIR/qc_imgs2 ]
then
    mkdir $SUBJECTS_DIR/qc_imgs2
fi

t1_struct=($(ls $SUBJECTS_DIR/T1_std.nii*))  # Get T1w data
#------------------------------------------------------


##############  FreeSurfer Segmentation  ##############
echo "bbb_pipe(): Segmenting brain structures using FastSurfer ..."



recon-all -subject freesurfer_output_on_T1_std -i $t1_struct -all 

mri_convert \
$SUBJECTS_DIR/freesurfer_output_on_T1_std/mri/aparc.DKTatlas+aseg.mgz \
$SUBJECTS_DIR/freesurfer_output_on_T1_std/mri/aparc+aseg.nii.gz

# Reslice FS segmentation to original T1 dims...
mri_convert \
-rl $t1_struct \
-rt nearest $SUBJECTS_DIR/freesurfer_output_on_T1_std/mri/aparc+aseg.nii.gz \
$SUBJECTS_DIR/freesurfer_output_on_T1_std/mri/aparc+aseg_reslice.nii.gz

# Restore ROIs from float to integer datatype
fslmaths \
$SUBJECTS_DIR/freesurfer_output_on_T1_std/mri/aparc+aseg_reslice.nii.gz \
$SUBJECTS_DIR/freesurfer_output_on_T1_std/mri/aparc+aseg_reslice_int32.nii.gz \
-odt int
	
echo "bbb_pipe(): PROCESSING COMPLETE"
}

# run the BBB_freesurfer_visualisation .py code to label and to see
# mri_convert $t1_struct \
# SUBJECTS_DIR/$subj_id/freesurfer_output_on_T1_std/mri/t1_struct.mgz --in_type nii --out_type mgz

# input=$SUBJECTS_DIR/$subj_id/freesurfer_output_on_T1_std/mri/t1_struct.mgz
# output=$SUBJECTS_DIR/$subj_id/freesurfer_output_on_T1_std/mri/aparc.DKTatlas+aseg.deep.mgz \

# freeview -v $input -v $output:colormap=lut


