# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import nibabel as nib

input_path = '/Users/Jayvik/Desktop'
subj = 'A22040401'


mask_file_path = 'T1_mask_binary/pred_mask_' + subj + '_T1.nii.gz'
T1_file_path = input_path+subj + '_T1.nii.gz'


mask = nib.load(mask_file_path)
mask_data = mask.get_fdata()

T1 = nib.load(T1_file_path)
T1_data = T1.get_fdata()

T1_masked_data = mask_data*T1_data

nifti = nib.Nifti1Image(T1_masked_data,T1.affine, T1.header)
nib.save(nifti, input_path + subj + 'T1_masked.nii.gz')