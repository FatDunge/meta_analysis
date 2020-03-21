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

    def caculate(self, studies, _model:str):
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

        if _model.lower() == 'fixed':
            m = model.FixedModel(effect_sizes, variances)
        elif _model.lower() == 'random':
            m = model.RandomModel(effect_sizes, variances)
        return m
