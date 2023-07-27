# Jayvik_fMRI

- Preprocessing Pipeline
  - AFNI
    - nifti_to_afni_conversion.txt: command line codes to convert between BRIK/HEAD and NIFTI files
    - template_afni_proc.py: afni's default pre-processing pipeline template structure, would probably work well with human brains
  - Mask_maker
    - model_2d_contrast_full_2.h5: deep learning model to generate masks for mice structural images
    - T1_RARE_mask_maker.py: code to run the masking model
    - binarizer.py: erode and dilate mask to clean it up
  - pipeline.py: full pre-processing pipeline from start to finish ***
  - random
    - collect_filenames.py: collects all filenames in a folder
    - graymatter.py: generates gray matter mask from csf, wm, and atlas masks
    - restack_slices.py: takes fMRI 3D Volumes and restacks into 4D NIFTI
    - separate_slices.py: takes fMRI 4D NIFTI and breaks down into 3D volumes
    - T1_resampler.py: resamples T1 images to 0.3mm^3 resolution

- Analysis
  - melodic.py: code to run FSL Melodic from terminal (important for mice since it can bypass masking step)
  - gift_setup.py: sets up files and filenames for Group ICA of fMRI Toolbox to organize better (important because subjects need to be in their own folder)
  - conn_setup.py: sets up files and filenames for Conn Toolbox (*** INCOMPLETE ***)
  - batch_test.m: runs Conn **
- NyquistGhost: document containing nyquist ghost research, papers, links to githubs of codes. (Minnesota code is of interest)
- Archive: unnecessary or incomplete scripts
  - conn_setup.py: sets up files and filenames for Conn Toolbox (*** INCOMPLETE ***)
  - apply_mask.py:   applies mask by multiplication between 2 nifti files (*** already part of pipeline.py ***)
  - echo_extracting.py:    extracts the first of 3 echos of a 4D image and removes the first 4 of those echos (*** already part of pipeline.py ***)
  - register_volumes_and_combine.py:   registers volumes using ANTs (*** already part of pipeline.py ***)
