#%%
import os
import nibabel as nib
import numpy as np
from meta_analysis import data, meta
def test_voxelmeta():
    path = '../test_data'
    mask = '../grey_matter_smoothed_005.nii'
    nii = nib.load(mask)
    mask = np.asarray(nii.dataobj)
    center_names = os.listdir(path)
    centers = []
    for center_name in center_names:
        center_path = os.path.join(path, center_name)
        files = os.listdir(center_path)
        group1 = []
        group2 = []
        for f in files:
            filepath = os.path.join(center_path, f)
            if '1_' in f:
                group1.append(filepath)
            elif '3_' in f:
                group2.append(filepath)
        g1 = data.ArrayGroup(label=1, datapathes=group1)
        g2 = data.ArrayGroup(label=3, datapathes=group2)
        center = data.Center(center_name, [g1, g2])
        centers.append(center)
    
    method = 'cohen_d'
    model = 'random'
    m = meta.VoxelMeta(centers)
    results = m.main(1,3,model,method,mask)

results = test_voxelmeta()
# %%
