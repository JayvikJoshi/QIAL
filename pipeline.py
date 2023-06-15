import nibabel as nib
import numpy as np
from scipy.ndimage import morphology
from nibabel import load, save, Nifti1Image, squeeze_image
import os
import sys, string, os
import pandas as pd


subj = 'A22040401'

my_path = '/Users/Jayvik/Desktop/preproc/Pipeline_Testing/'
input_path = my_path + subj + '/'

fmri_file_path = input_path + subj + '_fMRI.nii.gz'
bold = nib.load(fmri_file_path)
bold_data = bold.get_fdata()

T1_file_path = input_path + subj + '_T1.nii.gz'
T1 = nib.load(T1_file_path)
T1_data = T1.get_fdata()

mask_file_path = input_path + subj + '_mask.nii.gz'
mask = nib.load(mask_file_path)
mask_data = mask.get_fdata()

atlas_path = my_path + 'atlas.nii.gz'

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
os.system(f"3dTshift -tzero 0 -Fourier -prefix {STC_path} {func_echo1_path}")

VolReg_path = output + 'd_VolReg.nii.gz'
motion_measures_path = output + 'd_motion_measures.1D'
motion_transformation_matrix_path = output + 'd_motion_transformation_matrix.1D'
os.system(f"3dvolreg -verbose -zpad 1 -cubic -base 2 -1Dfile {motion_measures_path} -1Dmatrix_save {motion_transformation_matrix_path} -prefix {VolReg_path} {STC_path}")

T1_reoriented_path = output + 'e_T1_reoriented.nii.gz'
T1_to_Atlas_filename = 'e_T1_to_Atlas_'
T1_to_Atlas_matrix_path = output + f'{T1_to_Atlas_filename}0GenericAffine.mat'
T1_registered_path = output + 'e_T1_registered.nii.gz'
T1_resampled_path = output + 'e_T1_resampled.nii.gz'

os.system(f"c3d {T1_masked_path} -orient RAI -o {T1_reoriented_path}")
os.system(f"/opt/ANTs/bin/antsRegistration -v 1 -d 3 -m Mattes[ {atlas_path}, {T1_reoriented_path}, 1, 32, None ] -r [ {atlas_path}, {T1_reoriented_path}, 1 ]  -t affine[0.1] -c [300x300x0x0, 1e-8, 20 ] -s 4x2x1x0.5vox -f 6x4x2x1 -u 1 -z 1 -o {output}{T1_to_Atlas_filename}")
os.system(f"/opt/ANTs/bin/antsApplyTransforms -d 3 -e 0 --float -i {T1_reoriented_path} -r {atlas_path} -o {T1_registered_path} -t {T1_to_Atlas_matrix_path}")
os.system(f"/opt/ANTs/bin/ResampleImageBySpacing 3 {T1_registered_path} {T1_resampled_path} 0.1 0.1 0.1 0 0 0")

split_fMRI_path = output + 'split_fMRI/'
if not os.path.isdir(split_fMRI_path) : os.mkdir(split_fMRI_path)

for i in range(int(bold_data.shape[3])):
	bold_data_3d = bold_data[:,:,:,i]
	squeezed = squeeze_image(nib.Nifti1Image(bold_data_3d, bold.affine))
	nib.save(squeezed, split_fMRI_path + "vol_" + str(i) + ".nii.gz")

vol_0_path = split_fMRI_path + 'vol_0.nii.gz'
fMRI_to_T1_filename = 'f_fMRI_to_T1_'
fMRI_to_T1_matrix_path = split_fMRI_path + f'{fMRI_to_T1_filename}0GenericAffine.mat'

os.system("/opt/ANTs/bin/antsRegistration -v 1 -d 3 -m Mattes[ {T1_masked_path}, {vol_0_path}, 1, 32, None ] -r [ {T1_masked_path}, {vol_0_path}, 1 ]  -t affine[0.1] -c [300x300x0x0, 1e-8, 20 ] -s 4x2x1x0.5vox -f 6x4x2x1 -u 1 -z 1 -o {fMRI_to_T1_filename}")

fMRI_to_T1_registered_path = split_fMRI_path + f'{fMRI_to_T1_filename}registered_' #+(i).nii.gz

fMRI_to_Atlas_filename = 'f_fMRI_to_Atlas_'
fMRI_to_Atlas_reoriented_path = split_fMRI_path + f'{fMRI_to_Atlas_filename}reoriented_'
fMRI_to_Atlas_registered_path = split_fMRI_path + f'{fMRI_to_Atlas_filename}registered_'
fMRI_to_Atlas_resampled_path = split_fMRI_path + f'{fMRI_to_Atlas_filename}resampled_'

splitFile_suffix = ".nii.gz"

for i in range(int(bold_data.shape[3])):
	os.system(f"/opt/ANTs/bin/antsApplyTransforms -d 3 -e 0 --float -i {split_fMRI_path}vol_"+str(i)+f"{splitFile_suffix} -r {T1_masked_path} -o {fMRI_to_T1_registered_path}"+str(i)+f"{splitFile_suffix} -t {fMRI_to_T1_matrix_path}")
	os.system(f"c3d {fMRI_to_T1_registered_path}"+str(i)+f"{splitFile_suffix} -orient RAI -o {fMRI_to_Atlas_reoriented_path}"+str(i)+f"{splitFile_suffix}")
	os.system(f"/opt/ANTs/bin/antsApplyTransforms -d 3 -e 0 --float -i {fMRI_to_Atlas_reoriented_path}"+str(i)+f"{splitFile_suffix} -r {atlas_path} -o {fMRI_to_Atlas_registered_path}"+str(i)+f"{splitFile_suffix} -t {T1_to_Atlas_matrix_path}")
	os.system(f"/opt/ANTs/bin/ResampleImageBySpacing 3 {fMRI_to_Atlas_registered_path}"+str(i)+f"{splitFile_suffix} {fMRI_to_Atlas_resampled_path}"+str(i)+f"{splitFile_suffix} 0.3 0.3 0.3 0 0 0")

fMRI_to_T1_registered_recombined_path = output + 'f_fMRI_to_T1_registered_recombined.nii.gz'
fMRI_to_Atlas_reoriented_recombined_path = output + 'f_fMRI_to_Atlas_reoriented_recombined.nii.gz'
fMRI_to_Atlas_registered_recombined_path = output + 'f_fMRI_to_Atlas_registered_recombined.nii.gz'
fMRI_to_Atlas_resampled_recombined_path = output + 'f_fMRI_to_Atlas_resampled_recombined.nii.gz'

os.system(f'/opt/ANTs/bin/ImageMath 4 {fMRI_to_T1_registered_recombined_path} TimeSeriesAssemble 1 0 {fMRI_to_T1_registered_path}*{splitFile_suffix}')
os.system(f'/opt/ANTs/bin/ImageMath 4 {fMRI_to_Atlas_reoriented_recombined_path} TimeSeriesAssemble 1 0 {fMRI_to_Atlas_reoriented_path}*{splitFile_suffix}')
os.system(f'/opt/ANTs/bin/ImageMath 4 {fMRI_to_Atlas_registered_recombined_path } TimeSeriesAssemble 1 0 {fMRI_to_Atlas_registered_path}*{splitFile_suffix}')
os.system(f'/opt/ANTs/bin/ImageMath 4 {fMRI_to_Atlas_resampled_recombined_path } TimeSeriesAssemble 1 0 {fMRI_to_Atlas_resampled_path}*{splitFile_suffix}')

