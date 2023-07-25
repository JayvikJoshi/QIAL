import nibabel as nib
import numpy as np
from scipy.ndimage import morphology
from nibabel import load, save, Nifti1Image, processing
import os

input_path = '/Users/Jayvik/Desktop/smoothing_inputs/'
filepath='all_ica.nii'
file = nib.load(input_path + filepath)
file_data=file.get_fdata()
output = "/Users/Jayvik/Desktop/smoothed/"

files_list=os.listdir(input_path)   
for file in files_list:
    single_file=nib.load(input_path + file)
    single_file = nib.processing.smooth_image(single_file, (.25,.25,.25), mode='constant')
    file_data=single_file.get_fdata()
    #file_data[ file_data > 0 ] = 1
    file_result= nib.Nifti1Image(file_data, single_file.affine, single_file.header)
    nib.save(single_file, output+file[:-4]+"_smoothed.nii")
    