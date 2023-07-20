import nibabel as nib

# Load the masks
in_dir = "/Users/Jayvik/Desktop/"
atlas_path = "atlas_mask.nii"
csf_path = "csf_mask.nii"
wm_path = "wm_mask.nii"

atlas = nib.load(in_dir + atlas_path)
csf = nib.load(in_dir + csf_path)
wm = nib.load(in_dir + wm_path)

# Get the data from the masks
atlas_data = atlas.get_data()
csf_data = csf.get_data()
wm_data = wm.get_data()

# Subtract the masks
gm_data = atlas_data - csf_data - wm_data

# Save the result
nib.save(nib.Nifti1Image(gm_data, atlas.affine, atlas.header), in_dir + "gm_mask.nii")
