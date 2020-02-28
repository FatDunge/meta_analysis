from dataclasses import dataclass
from typing import Any

import nibabel as nib
import numpy as np

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

def cal_mean_std_n(arrays):
    mean = np.mean(arrays, axis=0)
    std = np.std(arrays, axis=0)
    n = arrays.shape[0]
    return mean, std, n

@dataclass
class Group(object):
    mean: Any
    std: Any
    count: int
    label: Any

    def __post_init__(self):
        self.mean = np.asarray(self.mean)
        self.std = np.asarray(self.std)

    def get_label(self):
        return self.label
    
    def get_mean_std_count(self, index):
        if index:
            return self.mean[index], self.std[index], self.count
        else:
            return self.mean.item(0), self.std.item(0), self.count

@dataclass
class Study(object):
    name: str
    method: str
    group_experimental: Group
    group_control: Group

    def __post_init__(self):
        self.s = None
        if self.method == 'cohen_d':
            self.func = self.cohen_d

    def cohen_d(self, index):
        m1, s1, n1 = self.group_experimental.get_mean_std_count(index)
        m2, s2, n2 = self.group_control.get_mean_std_count(index)
        self.s = np.sqrt(((n1-1)*(s1**2)+(n2-1)*(s2**2))/(n1+n2-2))
        d = (m1 - m2) / s
        return d

    def get_effect_size(self, index=None):
        return self.func(index)
    
    def get_variance(self):
        if not self.s:
            self.func()
        return self.s

class Center(object):
    def __init__(self, name, groups):
        self.name = name
        self.build_dic(groups)

    def build_dic(self, groups):
        group_dict = {}
        for group in groups:
            label = group.get_label()
            if label not in group_dict:
                group_dict[label] = group
            else:
                raise ValueError("Two Groups have same label, please check or merge.")
        self.group_dict = group_dict

    def gen_study(self, label1, label2, method):
        if label1 in self.group_dict and label2 in self.group_dict:
            return Study(self.name, method,
                         self.group_dict[label1],
                         self.group_dict[label2])

