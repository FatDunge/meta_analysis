""" data module mainly use to manage core data

Function:
    load_array(path): load nii file by filepath, return dataobj array
    load_arrays(pathes): load list of nii filepath, return ndarray of array
    cal_mean_std_n(arrays, axis): calculate arrays mean, std, n by axis

Class:
    Group(object): A specific label group within a center
    Study(object): A study which hold two group's mean, std, count,
                   can caculate its effect size, variance.
    Center(object): A center holds lots of group, generate study.

Author: Kang Xiaopeng
Data: 2020/02/20
E-mail: kangxiaopeng2018@ia.ac.cn

This file is part of meta_analysis.

meta_analysis is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

meta_analysis is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with meta_analysis.  If not, see <https://www.gnu.org/licenses/>.
"""

from dataclasses import dataclass
from typing import Any

import nibabel as nib
import numpy as np

from . import utils

class Group(object):
    """ Group of specific label, holds origin data or mean, std, count.

    Attributes:
        label: int or string, this group's label.
        datas: original datas, 1d array, used to caculate mean, std, count.
        mean: float, mean of datas
        std: float, std of datas
        count: int, count of datas
        shape: default 1, used to check mean's dimension.
        not_msn: bool, whether one of (mean, std, count) is None

    Function:
        check_property(): check if not_msn and not datas
        is_same_shape(): check if other group has same shape with self
        get_label(): return label
        get_mean_std_count(): return mean, std, count
    """

    def __init__(self, label, datas=None,
                 mean=None, std=None, count=None):
        """ init, Note must input (mean, std, count) or datas
        Args:
            label: int or string, this group's label.
            datas: original datas, 1d array, used to caculate mean, std, count.
            mean: float, mean of datas
            std: float, std of datas
            count: float, count of datas
        """
        super().__init__()
        self.label = label
        self.datas = datas
        self.mean = mean
        self.std = std
        self.count = count
        self.shape = 1

    def check_property(self):
        self.not_msn = self.mean is None or self.std is None or not self.count
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
                self.mean, self.std, self.count = utils.cal_mean_std_n(np.asarray(self.datas))
                self.not_msn = False
        return self.mean, self.std, self.count

@dataclass
class Study(object):
    """ meta analysis study, used to caculate effect size and variance.

    Attributes:
        name: str, study name
        method: str, method to caculate effect size 
        mean1: float, mean of experimental group
        std1: float, std of experimental group
        count1: int, count of experimental group
        mean2: float, mean of control group
        std2: float, std of control group
        count2: int, count of control group
        s: float, pooled standard deviation

    Function:
        cohen_d(): cohen'd effect size
        get_effect_size(): use specific method to return effect size
        get_variance(): return variance
    """
    name: str
    method: str
    group_experimental: Group
    group_control: Group

    def __post_init__(self):
        self.s = None
        if self.method == 'cohen_d':
            self.func = self.cohen_d
        elif self.method == 'hedge_g':
            self.func = self.hedge_g

    def cohen_d(self):
        """ details in https://en.wikipedia.org/wiki/Effect_size
        """
        m1, s1, n1 = self.group_experimental.get_mean_std_count()
        m2, s2, n2 = self.group_control.get_mean_std_count()
        s = np.sqrt(((n1-1)*(s1**2)+(n2-1)*(s2**2))/(n1+n2-2))
        d = (m1 - m2) / s
        self.s = s
        return d

    def hedge_g(self):
        g = self.cohen_d()
        n1, n2 = self.count1, self.count2
        g_star = (1-3/(4*(n1+n2)-9)) * g
        return g_star

    def get_effect_size(self):
        return self.func()

    def get_variance(self):
        if not self.s:
            self.func()
        return self.s ** 2

class Center(object):
    """ meta analysis study, used to caculate effect size and variance.

    Attributes:
        name: str, center name
        shape: group's mean shape, used to check center consistency
        group_dict: dict of Group instance, {name1: group1, ...}

    Function:
        check_groups(): check groups consistency
        build_dict(groups): build dict for better index
        check_label(label): check whether has group of specific label in this center
        gen_study(label1, label2, method, index): generate Study instance
        gen_region_study(label1, label2, method, mask, region_label): generate Region Study
    """
    def __init__(self, name, groups):
        self.name = name
        self.check_groups(groups)
        self.build_dict(groups)

    def check_groups(self, groups):
        #check shape
        base_group = groups[0]
        self.shape = base_group.shape
        class_name = type(base_group).__name__
        for group in groups:
            #check Group class
            assert type(group).__name__ == class_name
            assert base_group.is_same_shape(group)

    def is_same_shape(self, center):
        return self.shape == center.shape

    def build_dict(self, groups):
        group_dict = {}
        for group in groups:
            label = group.get_label()
            if label not in group_dict:
                group_dict[label] = group
            else:
                raise ValueError("Two Groups have same label, please check or merge.")
        self.group_dict = group_dict

    def check_label(self, label):
        if label not in self.group_dict:
            print('Couln\'t found [label:{}] group in [center:{}]'.format(label, self.name))
            return False
        return True

    def gen_study(self, label1, label2, method):
        if not self.check_label(label1):
            return
        elif not self.check_label(label2):
            return
        else:
            group_experimental = self.group_dict[label1]
            group_control = self.group_dict[label2]
            return Study(self.name, method,
                         group_experimental, group_control)

class Centers():
    def __init__(self, center_list):
        super().__init__()
        self.center_list = center_list

    def gen_studies(self, label1, label2, method):
        studies = []
        for center in self.center_list:
            study = center.gen_study(label1, label2, method)
            studies.append(study)
        return studies
