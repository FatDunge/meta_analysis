import nibabel as nib
import numpy as np
from .study import Study

def load_array(path):
    nii = nib.load(path)
    array = np.asarray(nii.dataobj)
    array = np.nan_to_num(array)
    return array

def load_arrays(pathes):
    arrays = np.array([])
    if pathes:
        arrays = np.stack([load_array(path) for path in pathes], axis=0)
    return arrays

def get_mean_std_n(arrays):
    mean = np.mean(arrays, axis=0)
    std = np.std(arrays, axis=0)
    n = arrays.shape[0]
    return mean, std, n

class Center(object):
    def __init__(self, name, filepathes1, filepathes2):
        self.name = name
        arrays1 = load_arrays(filepathes1)
        arrays2 = load_arrays(filepathes2)
        self.array_shape = arrays1[0].shape

        self.m1, self.s1, self.n1 = get_mean_std_n(arrays1)
        self.m2, self.s2, self.n2 = get_mean_std_n(arrays2)

    def gen_study(self, index):
        if self.n1 and self.n2:
            return Study(self.name,
                         self.m1[index], self.s1[index], self.n1[index],
                         self.m2[index], self.s2[index], self.n2[index])
    
    