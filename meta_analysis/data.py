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

class Group(object):
    def __init__(self, label, datas=None,
                 mean=None, std=None, count=None):
        super().__init__()
        self.label = label
        self.datas = datas
        self.mean = mean
        self.std = std
        self.count = count
        self.shape = 1

    def check_property(self):
        self.not_msn = not self.mean or not self.std or not self.count
        if self.not_msn:
            if self.datas:
                return True
            elif not self.datas:
                raise AttributeError('Need input for datas or (mean, std, count)')
        else:
            return True
    
    def is_same_shape(self, group):
        return self.shape == group.shape

    def get_label(self):
        return self.label
    
    def get_mean_std_count(self):
        if self.check_property():
            if self.not_msn:
                self.mean, self.std, self.count = cal_mean_std_n(np.asarray(self.datas))
                self.not_msn = False
        return self.mean, self.std, self.count

class ArrayGroup(Group):
    def __init__(self, label, datas=None,
                 mean=None, std=None, count=None,
                 datapathes=None):
        super().__init__(label, datas, mean, std, count)
        self.datapathes = datapathes
        self.gen_shape()
        
    def check_property(self):
        self.not_msn = not self.mean or not self.std or not self.count
        if self.not_msn:
            if self.datas is None:
                if self.datapathes:
                    self.datas = load_arrays(self.datapathes)
                    return True
                elif not self.datapathes:
                    raise AttributeError('Need at least one input for datas\
                                        or (mean, std, count) or datapathes')
            else:
                return True
        else:
            return True

    def gen_mean_std_count(self):
        if self.check_property():
            if self.not_msn:
                self.mean, self.std, self.count = super().get_mean_std_count()

    def gen_shape(self):
        self.gen_mean_std_count()
        self.shape = self.mean.shape

    def get_mean_std_count(self, index):
        return self.mean[index], self.std[index], self.count

@dataclass
class Study(object):
    name: str
    method: str
    mean1: float
    std1: float
    count1: int
    mean2: float
    std2: float
    count2: int

    def __post_init__(self):
        self.s = None
        if self.method == 'cohen_d':
            self.func = self.cohen_d

    def cohen_d(self):
        m1, s1, n1 = self.mean1, self.std1, self.count1
        m2, s2, n2 = self.mean2, self.std2, self.count2
        s = np.sqrt(((n1-1)*(s1**2)+(n2-1)*(s2**2))/(n1+n2-2))
        d = (m1 - m2) / s
        self.s = s
        return d

    def get_effect_size(self):
        return self.func()

    def get_variance(self):
        if not self.s:
            self.func()
        return self.s ** 2

class Center(object):
    def __init__(self, name, groups):
        self.name = name
        self.groups = groups
        self.check_groups()
        self.build_dic()

    def check_groups(self):
        base_group = self.groups[0]
        self.shape = base_group.shape
        class_name = type(base_group).__name__
        for group in self.groups:
            assert type(group).__name__ == class_name
            assert base_group.is_same_shape(group)

    def is_same_shape(self, center):
        return self.shape == center.shape

    def build_dic(self):
        group_dict = {}
        for group in self.groups:
            label = group.get_label()
            if label not in group_dict:
                group_dict[label] = group
            else:
                raise ValueError("Two Groups have same label, please check or merge.")
        self.group_dict = group_dict

    def gen_study(self, label1, label2, method, index=None):
        if label1 not in self.group_dict:
            print('Couln\'t found [label:{}] group in [center:{}]'.format(label1, self.name))
        elif label2 not in self.group_dict:
            print('Couln\'t found [label:{}] group in [center:{}]'.format(label2, self.name))
        else:
            if isinstance(self.group_dict[label1], ArrayGroup):
                if not index:
                    raise ValueError('ArrayGroup needs index to generate Study')
                else:
                    m1, s1, n1 = self.group_dict[label1].get_mean_std_count(index)
                    m2, s2, n2 = self.group_dict[label2].get_mean_std_count(index)                    
            else:
                m1, s1, n1 = self.group_dict[label1].get_mean_std_count()
                m2, s2, n2 = self.group_dict[label2].get_mean_std_count()
            return Study(self.name, method, m1, s1, n1, m2, s2, n2)
        
