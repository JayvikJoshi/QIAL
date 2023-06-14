

import nibabel as nib
import numpy as np
from scipy.ndimage import morphology
from nibabel import load, save, Nifti1Image, squeeze_image
import os
import sys, string, os
import pandas as pd


os.system("/opt/ANTs/bin/antsRegistration -v 1 -d 3 -m Mattes[ 01c_anat_masked.nii.gz, 0.nii.gz, 1, 32, None ] -r [ 01c_anat_masked.nii.gz, 0.nii.gz, 1 ]  -t affine[0.1] -c [300x300x0x0, 1e-8, 20 ] -s 4x2x1x0.5vox -f 6x4x2x1 -u 1 -z 1 -o 04c_")


for i in range(3):
	os.system("/opt/ANTs/bin/antsApplyTransforms -d 3 -e 0 --float -i "+str(i)+".nii.gz -r 01c_anat_masked.nii.gz -o 04c_"+str(i)+"_registered.nii.gz -t 04c_0GenericAffine.mat")
	os.system("c3d 04c_"+str(i)+"_registered.nii.gz -orient RAI -o 04d_"+str(i)+"_RAI.nii.gz")
	os.system("/opt/ANTs/bin/antsApplyTransforms -d 3 -e 0 --float -i 04d_"+str(i)+"_RAI.nii.gz -r 00_atlas.nii.gz -o 04d_"+str(i)+"_registered.nii.gz -t 04a_0GenericAffine.mat ")
	os.system("/opt/ANTs/bin/ResampleImageBySpacing 3 04d_"+str(i)+"_registered.nii.gz out_resampled/04d_"+str(i)+"_resampled.nii.gz 0.3 0.3 0.3 0 0 0")

os.system(f'/opt/ANTs/bin/ImageMath 4 /Users/Jayvik/Desktop/preproc/QC_testing/Test2/squeezed/combined/_4D.nii.gz TimeSeriesAssemble 1 0 /Users/Jayvik/Desktop/preproc/QC_testing/Test2/squeezed/out_resampled/04d_*_resampled.nii.gz')
