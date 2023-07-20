
"""
Created on Fri Feb  3 16:37:17 2023

@author: ali
"""

import nibabel as nib
import numpy as np
from scipy.ndimage import morphology
from nibabel import load, save, Nifti1Image, squeeze_image
import os
import sys, string, os
import pandas as pd
import openpyxl

mypath= '/Users/Jayvik/Desktop/regs/'

#construct csf and wm mask
#
list_fmri_folders = os.listdir(mypath)
list_of_subjs_long = [i for i in list_fmri_folders if 'regs' in i]


#os.system('/Applications/ANTS/antsApplyTransforms -d 3 -e 0 --float  -u float -i ' +mypath + 'olfac_amyg_hippo.nii.gz -n NearestNeighbor -r '+label_path_res+" -o "+mypath + 'olfac_amyg_hippo_0P3.nii.gz')

for subj in list_of_subjs_long:
    #print(subj)
    #fmri_file = list_fmir_folders_path +subj + "/ses-1/func/" + subj +"_ses-1_bold.nii.gz"
    #nib.load(fmri_file)
    #python_command = "python /mnt/munin6/Badea/Lab/mouse/fmri_pipeline/fmri_prep.py "+subj
    #job_name = job_descrp + "_"+ subj
    #command = GD + "submit_sge_cluster_job.bash " + sbatch_folder_path + " "+ job_name + " 0 0 '"+ python_command+"'"
    #os.system(command)
    subj_path_T1 = mypath + subj
    ref = mypath + '/chass_atlas_labels.nii'
    T1_res = mypath + subj[:9] + "_T1.nii.gz"
    os.system(f"/Users/alex/ANTS_install/build/ANTS-build/Examples/antsApplyTransforms -d 3 -e 0 --float  -u float -i {subj_path_T1} -n NearestNeighbor -r {ref} -o  {T1_res}")

'''

label_nii_wm_data =label_nii.get_fdata()*0

for wm in vol_index_wm:
    #print(csf)
    label_nii_wm_data[  data_label == int(wm)] = 1
    
    
    
file_result= nib.Nifti1Image(label_nii_wm_data, label_nii.affine, label_nii.header)
nib.save(file_result,mypath + 'chass_symmetric3/wm_mask.nii.gz'  )
os.system('/Applications/ANTS/antsApplyTransforms -d 3 -e 0 --float  -u float -i ' +mypath +'chass_symmetric3/wm_mask.nii.gz -n NearestNeighbor -r '+label_path_res+" -o "+mypath +'chass_symmetric3/wm_mask_0p3.nii.gz') 


########### making a mask out of labels

label_path= mypath +'/chass_symmetric3/chass_symmetric3_labels_PLI_res.nii.gz'
label_nii=nib.load(label_path)
mask_labels_data = label_nii.get_fdata()
mask_labels = np.unique(mask_labels_data)
mask_labels=np.delete(mask_labels, 0)
mask_of_label =label_nii.get_fdata()*0

for vol in mask_labels:
    mask_of_label[  mask_labels_data == int(vol)] = 1
    
file_result= nib.Nifti1Image(mask_of_label, label_nii.affine, label_nii.header)
nib.save(file_result,mypath + 'chass_symmetric3/mask_of_label.nii.gz'  )    
'''   
