import nibabel as nib
import os, sys, glob
import numpy as np

folder = '/Users/Jayvik/Desktop/nyquist/split3d/split2d_Nfix/'

files = glob.glob(os.path.join(folder,'*.nii'))
files = sorted(files)

filepath = '/Users/Jayvik/Desktop/nyquist/split3d/vol_0.nii.gz'
origstack_nii = nib.load(filepath)
newstack_array = np.zeros(origstack_nii.shape)

for i,file in enumerate(files):
	newstack_array[:,:,i] = nib.load(file).get_fdata()

newstack_nii = nib.Nifti1Image(newstack_array, origstack_nii.affine, origstack_nii.header)
nib.save(newstack_nii, os.path.join('/Users/Jayvik/Desktop/nyquist/split3d/stack3d_Nfix',os.path.basename(filepath)))
	