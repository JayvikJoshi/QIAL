import os

in_dir = r"C:\Users\jayvi\Desktop\QIAL\GIFT Testing\input\APOE4_4HN\\"

"""
for file in os.listdir(in_dir):
    filepath = in_dir + file
    if file.endswith(".nii.gz"):
        os.system(f"gzip.exe -d {filepath}")

    if file.endswith(".nii"):
        out_dir = in_dir + "APOE4_" + file[0:9] + "\\"
        new_filepath = out_dir + file
        if not os.path.isdir(out_dir) : os.mkdir(out_dir)
        os.rename(filepath, new_filepath)
"""


#for file in os.listdir(in_dir):
#    os.rename(in_dir + file, in_dir + "APOE4_" + file)    
