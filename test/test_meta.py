#%%
import os
import nibabel as nib
import numpy as np
from meta_analysis import main, mask, utils

def test_voxelmeta(_mask):
    path = r'../test_data'
    center_names = os.listdir(path)
    center_dict = {}
    for center_name in center_names:
        center_path = os.path.join(path, center_name)
        files = os.listdir(center_path)
        group1 = []
        group2 = []
        center = {}
        for f in files:
            filepath = os.path.join(center_path, f)
            if '1_' in f:
                group1.append(filepath)
            elif '3_' in f:
                group2.append(filepath)
        center[1] = group1
        center[3] = group2
        center_dict[center_name] = center

    method = 'cohen_d'
    model = 'random'
    results = main.voxelwise_meta_analysis(center_dict, 1, 3,
                                           model=model,method=method,mask=_mask)
    return results

mask_path = r'../grey_matter_smoothed_005.nii'
mask_nii = nib.load(mask_path)
_mask = np.asarray(mask_nii.dataobj)
_mask = mask.Mask(_mask)
results = test_voxelmeta(_mask)

result_names = ['es','var', 'se', 'll','ul','q','z','p']
output = r'../result/{}'
for result, name in zip(results, result_names):
    path = output.format(name)
    utils.gen_nii(result, mask_nii, path)
