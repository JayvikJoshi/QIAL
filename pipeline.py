import nibabel as nib
import numpy as np
from scipy.ndimage import morphology
from nibabel import load, save, Nifti1Image, squeeze_image
import os
import sys, string, os
import pandas as pd


subj = 'A22040401'

input_path = '/Users/Jayvik/Desktop/preproc/Pipeline_Testing/' + subj + '/'

fmri_file_path = input_path + subj + '_fMRI.nii.gz'
bold = nib.load(fmri_file_path)
bold_data = bold.get_fdata()

T1_file_path = input_path + subj + '_T1.nii.gz'
T1 = nib.load(T1_file_path)
T1_data = T1.get_fdata()

mask_file_path = input_path + subj + '_mask.nii.gz'
mask = nib.load(mask_file_path)
mask_data = mask.get_fdata()

output = input_path + 'output/'
if not os.path.isdir(output) : os.mkdir(output)

T1_masked_data = mask_data * T1_data

T1_masked_nii=nib.Nifti1Image(T1_masked_data, T1.affine, T1.header)
nib.save(T1_masked_nii, output + 'a_T1_masked.nii.gz')
T1_masked_path= output + 'a_T1_masked.nii.gz'

beginVolume = 4
seq = range(beginVolume*3, bold_data.shape[3], 3)
echo1 = bold_data[:,:,:,seq]

func_echo1_nii = nib.Nifti1Image(echo1,bold.affine)
nib.save(func_echo1_nii, output + 'b_func_echo1.nii.gz')
func_echo1_path = output + 'b_func_echo1.nii.gz'

STC_path = output + 'c_STC.nii.gz'
os.system(f'3dTshift -tzero 0 -Fourier -prefix {STC_path} {func_echo1_path}')

VolReg_path = output + 'd_VolReg.nii.gz'
motion_measures_path = output + 'd_motion_measures.1D'
motion_transformation_matrix_path = output + 'd_motion_transformation_matrix.1D'
os.system(f'3dvolreg -verbose -zpad 1 -cubic -base 2 -1Dfile {motion_measures_path} -1Dmatrix_save {motion_transformation_matrix_path} -prefix {VolReg_path} {STC_path}')

