""" meta module as meta analysis entry point

Class:
    Meta(object): perform meta analysis
    VoxelMeta(Meta): perform voxel wise meta analysis

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
import nibabel as nib
import numpy as np

from . import model
from . import mask

class Meta(object):
    """Basic class to perform meta analysis

    Function:
        caculate(studies, model): caculate effect size and variance of each study,
                                  feed them into model, return model
    """
    def __init__(self):
        super().__init__()

    def caculate(self, studies, model:str):
        """caculate effect size and variance, feed to model
        Args:
            studies: list of Study
            model: 'fixed' or 'random'
        Return:
            m: instance of Model
        """
        effect_sizes = []
        variances = []
        for study in studies:
            es = study.get_effect_size()
            v = study.get_variance()

            effect_sizes.append(es)
            variances.append(v)

        if model.lower() == 'fixed':
            m = model.FixedModel(effect_sizes, variances)
        elif model.lower() == 'random':
            m = model.RandomModel(effect_sizes, variances)
        return m

class VoxelMeta(Meta):
    """perform voxelwise meta analysis

    Attributes:
        centers: list of Center instance
        results_array: result's array

    Function:
        check_centers(): check centers consistency
        main(): generate studies and caculate results per voxel
    """
    def __init__(self, centers):
        super().__init__()
        self.centers = centers
        self.check_centers()
        self.results_array = None

    def check_centers(self):
        """check whether centers datas has same shape
        """
        base_center = self.centers[0]
        for center in self.centers:
            assert base_center.is_same_shape(center)

    def main(self, label1, label2, model:str, method:str, _mask=None):
        """ for every voxel generate studies the caculate
        Args:
            label1: experimental group's label
            label2: control group's label
            model: 'fixed' or 'random'
            method: effect size method. 'cohen_d'
            _mask: Mask instance, None means all voxel.
        Return:
            results_array: ndarray, concatenate results return from model,
                           shape=(8, array_shape),
                           [0]: total effect size
                           [1]: total variance
                           [2]: standard error
                           [3]: 95% lower limits
                           [4]: 95% upper limits
                           [5]: heterogeneity
                           [6]: z test value
                           [7]: p value of z value
        """
        array_shape = self.centers[0].shape
        if _mask is None:
            array = np.ones(array_shape)
            _mask = mask.Mask(array)
        assert _mask.get_shape() == array_shape
        i = 0
        # FIXME: time issue(each loop takes 0.15s), try split code into gen_study and caculate
        nonzero_index = _mask.get_nonzero_index()
        for index in nonzero_index:
            index = tuple(index)
            # construct studies
            studies = []
            for center in self.centers:
                study = center.gen_study(label1, label2, method, index)
                studies.append(study)
            m = super().caculate(studies, model)
            results = m.get_results()

            if i == 0:
                results_array = np.zeros(((len(results),) + array_shape))
                i += 1
            for j in range(len(results)):
                results_array[j][index] = results[j]
        self.results_array = results_array
        return self.results_array

