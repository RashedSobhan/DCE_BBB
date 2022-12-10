#!/bin/bash
bbb_pipe(){
#=========================================================================
# TITLE: "BBB Segementation Pipeline, Version 0.1" (bbb_fastsurfer.sh)
# PURPOSE: Script for obtaining ROI from which signal enhancement will be extracted. 
# this code Performs segmentation of T1-weighted images in FreeSurfer. 
#
# INPUT: 
# AUTHORS: Donnie Cameron
# CREATED: 2022/05/05
# LAST UPDATED:  
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
subj_id="$@"



#source activate RS_conda_fastsurfer
# Source Free- and FastSurfer to allow function to be called here:
#export FREESURFER_HOME=/gpfs/software/ada/freesurfer/6.0

#export FREESURFER_HOME=/gpfs/home/wae16sru/freesurfer
#source $FREESURFER_HOME/SetUpFreeSurfer.sh
#export PATH=$PATH:/gpfs/software/ada/freesurfer/6.0/mni/bin
#export PERL5LIB=$PERL5LIB:/gpfs/home/wae16sru/mni-perllib

#module add freesurfer/6.0
#setupfreesurfer
#module load FastSurfer/0.0

module add fsl
#module load python/anaconda/2020.10/3.8
module add freesurfer
export SUBJECTS_DIR=~/BBB_sample_data/$subj_id
. /gpfs/software/ada/freesurfer/6.0/FreeSurferEnv.sh

module add FastSurfer


#echo "enter subject ID: "
#read x 
#dir= "/gpfs/home/wae16sru/BBB_sample_data/$x" 

echo "bbb_pipe(): Dataset (\$@) passed - \"$@\""



t1_struct=($(ls $SUBJECTS_DIR/T1_std.nii.gz))  # Get T1w data
#------------------------------------------------------

# Create directories for FSL products and QC images
if [ ! -d $SUBJECTS_DIR/products2 ]
then
    mkdir $SUBJECTS_DIR/products2
fi
if [ ! -d $SUBJECTS_DIR/qc_imgs2 ]
then
    mkdir $SUBJECTS_DIR/qc_imgs2
fi

##############  FastSurfer Segmentation  ##############
echo "bbb_pipe(): Segmenting brain structures using FastSurfer ..."

cd /gpfs/software/ada/FastSufer/0.0/FastSurfer-master/  
./run_fastsurfer.sh \
--t1 $t1_struct \
--sid fastsurfer_output_on_T1_std \
--sd $SUBJECTS_DIR \
--seg_only \
--parallel \
--threads 8 \
--py python3.8

mri_convert \
$SUBJECTS_DIR/fastsurfer_output_on_T1_std/mri/aparc.DKTatlas+aseg.deep.mgz \
$SUBJECTS_DIR/fastsurfer_output_on_T1_std/mri/aparc+aseg.nii.gz

# Reslice FS segmentation to original T1 dims...
mri_convert \
-rl $t1_struct \
-rt nearest $SUBJECTS_DIR/fastsurfer_output_on_T1_std/mri/aparc+aseg.nii.gz \
$SUBJECTS_DIR/fastsurfer_output_on_T1_std/mri/aparc+aseg_reslice.nii.gz

# Restore ROIs from float to integer datatype
fslmaths \
$SUBJECTS_DIR/fastsurfer_output_on_T1_std/mri/aparc+aseg_reslice.nii.gz \
$SUBJECTS_DIR/fastsurfer_output_on_T1_std/mri/aparc+aseg_reslice_int32.nii.gz \
-odt int
	
echo "bbb_pipe(): PROCESSING COMPLETE"
}

# run the BBB_freesurfer_visualisation .py code to label and to see
# mri_convert $t1_struct \
# SUBJECTS_DIR/$subj_id/freesurfer_output/mri/t1_struct.mgz --in_type nii --out_type mgz

# input=$SUBJECTS_DIR/$subj_id/freesurfer_output/mri/t1_struct.mgz
# output=$SUBJECTS_DIR/$subj_id/freesurfer_output/mri/aparc.DKTatlas+aseg.deep.mgz \

# freeview -v $input -v $output:colormap=lut


