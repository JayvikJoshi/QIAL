import os

input_dir = '/Users/Jayvik/Desktop/run01_inputs/'

T1_dir = input_dir + 'processed_T1/'
fMRI_dir = input_dir + 'processed_fMRI/'

listdir_T1 = os.listdir(T1_dir)
listdir_fMRI = os.listdir(fMRI_dir)

list_T1 = [i for i in listdir_T1]
list_fMRI = [i for i in listdir_T1]

for subj in list_T1:
    subj_renamed =
