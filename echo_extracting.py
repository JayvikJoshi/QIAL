# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import nibabel as nib

beginVolume = 4

subj = 'A22040401'
input_path = '/Users/Jayvik/Desktop/'
fmri_file_path=input_path+subj+'_fMRI.nii.gz'

bold=nib.load(fmri_file_path)
bold_data=bold.get_fdata()

seq = range(start*3, bold_data.shape[3], 3)
echo_1 = bold_data[:,:,:,seq]


nifti = nib.Nifti1Image(echo1,bold.affine)
nib.save(nifti, input_path + subj + 'echo1_.nii.gz')

