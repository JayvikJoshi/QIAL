
import nibabel as nib
from nibabel import load, save, Nifti1Image, squeeze_image
import numpy as np
import pandas as pd
import sys, string, os


fmri_file_path = '/Users/Jayvik/Desktop/preproc/QC_testing/Test2/03b_func_volreg.nii.gz'
bold = nib.load(fmri_file_path)
bold_data = bold.get_fdata()

squeezed_path = '/Users/Jayvik/Desktop/preproc/QC_testing/Test2/squeezed'
if not os.path.isdir(squeezed_path) : os.mkdir(squeezed_path)

for i in range(int(bold_data.shape[3])):
	bold_data_3d = bold_data[:,:,:,i]
	squeezed=squeeze_image(nib.Nifti1Image(bold_data_3d, bold.affine))
	nib.save(squeezed, "/Users/Jayvik/Desktop/preproc/QC_testing/Test2/squeezed/" + str(i) + ".nii.gz")
	
