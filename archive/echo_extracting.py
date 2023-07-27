# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import nibabel as nib

beginVolume = 4

subj = ''
input_path = '/Users/Jayvik/Desktop/preproc/QC_testing/Test2/'
fmri_file_path=input_path+subj+'00_func.nii.gz'

bold=nib.load(fmri_file_path)
bold_data=bold.get_fdata()

seq = range(beginVolume*3, bold_data.shape[3], 3)
echo_1 = bold_data[:,:,:,seq]


nifti = nib.Nifti1Image(echo_1,bold.affine)
nib.save(nifti, input_path + subj + '02_func_echo1.nii.gz')

