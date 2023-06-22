import nibabel as nib
import numpy as np
from scipy.ndimage import morphology
from nibabel import load, save, Nifti1Image, squeeze_image
import os
import sys, string, os
import pandas as pd

###inputs

subj = 'A22040401'
my_path = '/Users/Jayvik/Desktop/preproc/Pipeline_Testing/'
input_path = my_path + subj + '/'

afni_path = '/Users/Jayvik/abin/'
fsl_path = '/usr/local/fsl/share/fsl/bin/'
ants_path = '/opt/ANTs/bin/'
c3d_path = '/Applications/Convert3DGUI.app/Contents/bin/'

#inputting files
fMRI_file_path = input_path + subj + '_fMRI.nii.gz'
T1_file_path = input_path + subj + '_T1.nii.gz'
B0_file_path = input_path + subj + '_B0.nii.gz'
mask_file_path = input_path + subj + '_mask.nii.gz'
atlas_path = my_path + 'atlas.nii.gz'
atlas_mask_path = my_path + 'atlas_mask.nii.gz'
csf_mask_path = my_path + 'csf_mask.nii.gz'
wm_mask_path = my_path + 'wm_mask.nii.gz'





#loading fMRI, T1, B0, binarized mask data
fMRI = nib.load(fMRI_file_path)
fMRI_data = fMRI.get_fdata()
T1 = nib.load(T1_file_path)
T1_data = T1.get_fdata()
B0 = nib.load(B0_file_path)
B0_data = B0.get_fdata()
mask = nib.load(mask_file_path)
mask_data = mask.get_fdata()


#creating processing block folders

pb01 = input_path + 'pb01_T1_masking/'
if not os.path.isdir(pb01) : os.mkdir(pb01)

pb02 = input_path + 'pb02_Echo_Extraction/'
if not os.path.isdir(pb02) : os.mkdir(pb02)

pb03 = input_path + 'pb03_bias_correction/'
if not os.path.isdir(pb03) : os.mkdir(pb03)

pb04 = input_path + 'pb04_STC_VolReg/'
if not os.path.isdir(pb04) : os.mkdir(pb04)

pb06 = input_path + 'pb06_T1_Atlas_Reg/'
if not os.path.isdir(pb06) : os.mkdir(pb06)

pb07 = input_path + 'pb07_Vol_Atlas_Reg/'
if not os.path.isdir(pb07) : os.mkdir(pb07)

pb08 = input_path + 'pb08_Despike_Detrend_Blur_Scale/'
if not os.path.isdir(pb08) : os.mkdir(pb08)

pb09 = input_path + 'pb09_Denoise_Bandpass/'
if not os.path.isdir(pb09) : os.mkdir(pb09)


#masking T1 data and saving

T1_masked_data = mask_data * T1_data
T1_masked_nii=nib.Nifti1Image(T1_masked_data, T1.affine, T1.header)
nib.save(T1_masked_nii, pb01 + 'T1_masked.nii.gz')
T1_masked_path= pb01 + 'T1_masked.nii.gz'

#extracting echos and removing first 3 volumes

beginVolume = 4
seq = range(beginVolume*3, fMRI_data.shape[3], 3)
echo1 = fMRI_data[:,:,:,seq]
func_echo1_nii = nib.Nifti1Image(echo1,fMRI.affine)
nib.save(func_echo1_nii, pb02 + 'func_echo1.nii.gz')
func_echo1_path = pb02 + 'func_echo1.nii.gz'

#bias correction

#slice time correction and volume registration

STC_path = pb04 + 'STC.nii.gz'
os.system(f"{afni_path}3dTshift -tzero 0 -Fourier -prefix {STC_path} {func_echo1_path}")

VolReg_path = pb04 + 'VolReg.nii.gz'
motion_measures_path = pb04 + 'motion_measures.1D'
motion_transformation_matrix_path = pb04 + 'motion_transformation_matrix.1D'
os.system(f"{afni_path}3dvolreg -verbose -zpad 1 -cubic -base 2 -1Dfile {motion_measures_path} -1Dmatrix_save {motion_transformation_matrix_path} -prefix {VolReg_path} {STC_path}")

#T1 registering to Atlas

T1_reoriented_path = pb06 + 'T1_reoriented.nii.gz'
os.system(f"{c3d_path}c3d {T1_masked_path} -orient RAI -o {T1_reoriented_path}")

T1_to_Atlas_filename = 'T1_to_Atlas_'
T1_to_Atlas_matrix_path = pb06 + f'{T1_to_Atlas_filename}0GenericAffine.mat'
os.system(f"{ants_path}antsRegistration -v 1 -d 3 -m Mattes[ {atlas_path}, {T1_reoriented_path}, 1, 32, None ] -r [ {atlas_path}, {T1_reoriented_path}, 1 ]  -t affine[0.1] -c [300x300x0x0, 1e-8, 20 ] -s 4x2x1x0.5vox -f 6x4x2x1 -u 1 -z 1 -o {pb06}{T1_to_Atlas_filename}")

T1_registered_path = pb06 + 'T1_registered.nii.gz'
os.system(f"{ants_path}antsApplyTransforms -d 3 -e 0 --float -i {T1_reoriented_path} -r {atlas_path} -o {T1_registered_path} -t {T1_to_Atlas_matrix_path}")

T1_resampled_path = pb06 + 'T1_resampled.nii.gz'
os.system(f"{ants_path}ResampleImageBySpacing 3 {T1_registered_path} {T1_resampled_path} 0.1 0.1 0.1 0 0 0")

#splitting fMRI volumes

split_fMRI_path = pb07 + 'split_fMRI/'
if not os.path.isdir(split_fMRI_path) : os.mkdir(split_fMRI_path)

split = nib.load(VolReg_path)
split_data = split.get_fdata()

for i in range(int(split_data.shape[3])):
    split_data_i = split_data[:,:,:,i]
    squeezed = squeeze_image(nib.Nifti1Image(split_data_i, split.affine))
    nib.save(squeezed, split_fMRI_path + "vol_" + str(i) + ".nii.gz")

#fMRI registering to Atlas

vol_0_path = split_fMRI_path + 'vol_0.nii.gz'

fMRI_to_T1_filename = 'fMRI_to_T1_'
fMRI_to_T1_matrix_path = pb07 + f'{fMRI_to_T1_filename}0GenericAffine.mat'
os.system(f"{ants_path}antsRegistration -v 1 -d 3 -m Mattes[ {T1_masked_path}, {vol_0_path}, 1, 32, None ] -r [ {T1_masked_path}, {vol_0_path}, 1 ]  -t affine[0.1] -c [300x300x0x0, 1e-8, 20 ] -s 4x2x1x0.5vox -f 6x4x2x1 -u 1 -z 1 -o {pb07}{fMRI_to_T1_filename}")

#b0 registered to first volume

B0_to_fMRI_filename = 'B0_to_fMRI_'
B0_to_fMRI_registered_0_path = pb07 + f'{B0_to_fMRI_filename}registered_0.nii.gz'
os.system(f"{ants_path}antsApplyTransforms -d 3 -e 0 -i {B0_file_path} -r {vol_0_path} -u float -o {B0_to_fMRI_registered_0_path}")

B0_to_fMRI_applied_path = split_fMRI_path + f'{B0_to_fMRI_filename}applied_'
N4Bias_applied_path = split_fMRI_path + 'N4Bias_applied_'
bias_bold_path = split_fMRI_path + 'bias_bold'
fMRI_to_T1_registered_path = split_fMRI_path + f'{fMRI_to_T1_filename}registered_' #+(i).nii.gz
fMRI_to_Atlas_filename = 'fMRI_to_Atlas_'
fMRI_to_Atlas_reoriented_path = split_fMRI_path + f'{fMRI_to_Atlas_filename}reoriented_'
fMRI_to_Atlas_registered_path = split_fMRI_path + f'{fMRI_to_Atlas_filename}registered_'
fMRI_to_Atlas_resampled_path = split_fMRI_path + f'{fMRI_to_Atlas_filename}resampled_'
splitFile_suffix = ".nii.gz"

for i in range(75): #split_data.shape[3]
    #registered b0 applied to all volumes
    os.system(f"{fsl_path}fugue -i {split_fMRI_path}vol_"+str(i)+f"{splitFile_suffix} --dwell=0.00139 --loadfmap={B0_to_fMRI_registered_0_path} -s 0.3 --unwarpdir=x- -v -u {B0_to_fMRI_applied_path}"+str(i)+f"{splitFile_suffix}")
    #N4BiasFieldCorrection
    os.system(f"{ants_path}N4BiasFieldCorrection -d 3 -v 1 -s 2 -b [ 2x2x2 , 3 ] -c [ 100 x 50, 1e-6 ] -i {B0_to_fMRI_applied_path}"+str(i)+f"{splitFile_suffix} -o [{N4Bias_applied_path}"+str(i)+f"{splitFile_suffix},{bias_bold_path}"+str(i)+f"{splitFile_suffix}]")
    #all volumes registered to T1
    os.system(f"{ants_path}antsApplyTransforms -d 3 -e 0 --float -i {N4Bias_applied_path}"+str(i)+f"{splitFile_suffix} -r {T1_masked_path} -o {fMRI_to_T1_registered_path}"+str(i)+f"{splitFile_suffix} -t {fMRI_to_T1_matrix_path}")
    #all volumes registered to atlas (using T1 to atlas transformations)
    os.system(f"{c3d_path}c3d {fMRI_to_T1_registered_path}"+str(i)+f"{splitFile_suffix} -orient RAI -o {fMRI_to_Atlas_reoriented_path}"+str(i)+f"{splitFile_suffix}")
    os.system(f"{ants_path}antsApplyTransforms -d 3 -e 0 --float -i {fMRI_to_Atlas_reoriented_path}"+str(i)+f"{splitFile_suffix} -r {atlas_path} -o {fMRI_to_Atlas_registered_path}"+str(i)+f"{splitFile_suffix} -t {T1_to_Atlas_matrix_path}")
    os.system(f"{ants_path}ResampleImageBySpacing 3 {fMRI_to_Atlas_registered_path}"+str(i)+f"{splitFile_suffix} {fMRI_to_Atlas_resampled_path}"+str(i)+f"{splitFile_suffix} 0.3 0.3 0.3 0 0 0")


#recombining and saving volumes
fMRI_bias_corrected_recombined_path = pb07 + 'fMRI_bias_corrected_recombined.nii.gz'
os.system(f"{ants_path}ImageMath 4 {fMRI_bias_corrected_recombined_path} TimeSeriesAssemble 1 0 {N4Bias_applied_path}*{splitFile_suffix}")

fMRI_to_T1_registered_recombined_path = pb07 + 'fMRI_to_T1_registered_recombined.nii.gz'
os.system(f"{ants_path}ImageMath 4 {fMRI_to_T1_registered_recombined_path} TimeSeriesAssemble 1 0 {fMRI_to_T1_registered_path}*{splitFile_suffix}")

fMRI_to_Atlas_reoriented_recombined_path = pb07 + 'fMRI_to_Atlas_reoriented_recombined.nii.gz'
os.system(f"{ants_path}ImageMath 4 {fMRI_to_Atlas_reoriented_recombined_path} TimeSeriesAssemble 1 0 {fMRI_to_Atlas_reoriented_path}*{splitFile_suffix}")

fMRI_to_Atlas_registered_recombined_path = pb07 + 'fMRI_to_Atlas_registered_recombined.nii.gz'
os.system(f"{ants_path}ImageMath 4 {fMRI_to_Atlas_registered_recombined_path} TimeSeriesAssemble 1 0 {fMRI_to_Atlas_registered_path}*{splitFile_suffix}")

fMRI_to_Atlas_resampled_recombined_path = pb07 + 'fMRI_to_Atlas_resampled_recombined.nii.gz'
os.system(f"{ants_path}ImageMath 4 {fMRI_to_Atlas_resampled_recombined_path} TimeSeriesAssemble 1 0 {fMRI_to_Atlas_resampled_path}*{splitFile_suffix}")


#despiking

fMRI_despiked_path = pb08 + 'fMRI_despiked.nii.gz'
os.system(f"{afni_path}3dDespike -NEW -nomask -prefix {fMRI_despiked_path} {fMRI_to_Atlas_resampled_recombined_path}")

#detrending

#fMRI_detrend_path = pb08 + 'fMRI_detrend.nii.gz'
#os.system(f"{afni_path}3dDetrend -polort 9 -prefix {fMRI_detrend_path} {fMRI_despiked_path}")

#blurring

fMRI_blur_path = pb08 + 'fMRI_blur.nii.gz'
os.system(f"{afni_path}3dBlurInMask -input {fMRI_despiked_path} -prefix {fMRI_blur_path} -FWHM 0.45 -mask {atlas_mask_path}")

#scaling

fMRI_3dTstat_mean_path = pb08 + 'fMRI_3dTstat_mean.nii.gz'
fMRI_scaled_path = pb08 + 'fMRI_scaled.nii.gz'
os.system(f"{afni_path}3dTstat -prefix {fMRI_3dTstat_mean_path} {fMRI_blur_path}")
os.system(f"{afni_path}3dcalc -a {fMRI_blur_path} -b {fMRI_3dTstat_mean_path} -c {atlas_mask_path} -expr 'c * min(200, a/b*100)*step(a)*step(b)' -prefix {fMRI_scaled_path} ")

#csf and wm denoising

csf_mean_path = pb09 + 'csf_mean.1D'
wm_mean_path = pb09 + 'wm_mean.1D'

os.system(f"{afni_path}3dTstat -mean -mask {csf_mask_path} -prefix {csf_mean_path} {fMRI_scaled_path}")
os.system(f"{afni_path}3dTstat -mean -mask {wm_mask_path} -prefix {wm_mean_path} {fMRI_scaled_path}")

X_xmat_1D = pb09 + 'X.xmat.1D'
X_jpg = pb09 +'X.jpg'
fitts_subj  = pb09 + 'fitts'
errts_subj =  pb09 + 'errts'
stats_subj = pb09 + 'stats'


os.system(f"{afni_path}3dDeconvolve -input {fMRI_scaled_path} -ortvec {csf_mean_path} {csf_mask_path} -ortvec {wm_mean_path} {wm_mask_path} -polort 5 -float -num_stimts 0 -fout -tout -x1D {X_xmat_1D} -xjpeg {X_jpg} -fitts {fitts_subj} -errts {errts_subj} -x1D_stop -bucket {stats_subj} -overwrite")

#bandpassing

errts_nii_gz = pb09 + subj + '_errts.nii.gz'
os.system(f"{afni_path}3dTproject -polort 0 -input {fMRI_scaled_path} -ort {X_xmat_1D} -passband 0.01 0.1 -prefix {errts_nii_gz}")
