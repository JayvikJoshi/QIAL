# Jayvik_fMRI

AFNI
- nifti_to_afni_conversion: code to convert from .nii.gz to BRIK/HEAD format and vice versa
- template_afni_proc.py: afni preprocessing template for resting state fMRI without using FreeSurfer

mask_maker
- T1_RARE_mask_maker.py: deep learning model for mask making
- binarizer.py: binarizes, erodes, and dilates mask to get rid of holes

Preprocessing Pipeline
- apply_mask.py: uses subject and mask input and multiples them
- echo_extracting.py: extracts first echo from fMRI dataset and then removes first 4 volumes
- register_volumes_and_combine.py: registers T1 to atlas, splits fMRI volumes, registers them to atlas as well, and then concatenates them
- separate_slices.py: separates volumes of fMRI data
- pipeline.py: runs processing on data (masking, echo extraction, slice time correction, volume registration, B0 warping, N4BiasFieldcorrection, registration to atlas, despiking, blurring, scaling, csf and wm denoising, bandpassing)

Analysis
- melodic.py: creates txt file with all path directories of .nii.gz errts files and runs FSL MELODIC, runs cluster function (by splitting each component output, running the function on each, and concatenating), and concatenates thresh_zstat output
