import nibabel as nib
import numpy as np
from scipy.ndimage import morphology
from nibabel import load, save, Nifti1Image, squeeze_image
import os
import sys, string, os
import pandas as pd
import shutil
###inputs

test = 'temp'
my_path = '/Users/Jayvik/Desktop/analysis/FSL_MELODIC/'
input_path = my_path + test + '/'

fsl_path = '/usr/local/fsl/share/fsl/bin/'
ants_path = '/opt/ANTs/bin/'

errts_list = ''
errts_list_path = input_path +'errts_list'
fsl_output_path = input_path + 'output/'

final_path = input_path + 'final/'
if not os.path.isdir(final_path) : os.mkdir(final_path)

###RUN MELODIC

for filename in os.listdir(input_path):
    if filename.endswith(".nii.gz"):
        errts_list += input_path + filename + '\n'
    
with open (input_path+'errts_list', 'w') as text_file:
    text_file.write(errts_list)

melodic_path = fsl_output_path + 'melodic_IC.nii.gz'
os.system(f"{fsl_path}melodic -i {errts_list_path} -o {fsl_output_path} --nobet --mmthresh=0.5 --tr=1 --Oall")
shutil.copy(melodic_path, final_path + "melodic_out.nii.gz")

###RUN CLUSTER

split_melodic_path = fsl_output_path + 'split_melodic/'
if not os.path.isdir(split_melodic_path) : os.mkdir(split_melodic_path)

melodic = nib.load(melodic_path)
melodic_data = melodic.get_fdata()

for i in range(melodic_data.shape[3]):
    melodic_data_i = melodic_data[:,:,:,i]
    squeezed = squeeze_image(nib.Nifti1Image(melodic_data_i, melodic.affine))
    nib.save(squeezed, split_melodic_path + "melodic" + str(i) + ".nii.gz")
    os.system(f"{fsl_path}cluster -i {split_melodic_path}melodic" + str(i) + f".nii.gz -o {split_melodic_path}cluster" + str(i) + f".nii.gz --zthresh=1.645")
    
os.system(f"{ants_path}ImageMath 4 {final_path}cluster_out.nii.gz TimeSeriesAssemble 1 0 {split_melodic_path}cluster*.nii.gz")

###THRESH_ZSTAT + PROBMAP

stat_path = fsl_output_path + 'stats/'

os.system(f"{ants_path}ImageMath 4 {final_path}thresh_zstat_out.nii.gz TimeSeriesAssemble 1 0 {stat_path}thresh_zstat*.nii.gz")
os.system(f"{ants_path}ImageMath 4 {final_path}probmap_out.nii.gz TimeSeriesAssemble 1 0 {stat_path}probmap_*.nii.gz")
