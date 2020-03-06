""" Some helper function.

Function:
    gen_nii(array, template_nii, path): generate nii file using template's header and affine

Author: Kang Xiaopeng
Data: 2020/03/06
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
import os

import nibabel as nib

def gen_nii(array, template_nii, path=None):
    """ generate nii file using template's header and affine
        if input path then save nii in disk.
    Args:
        studies: list of Study
        model: 'fixed' or 'random'
    Return:
        m: instance of Model
    """
    affine = template_nii.affine
    header = template_nii.header
    nii = nib.Nifti1Image(array, affine, header)
    if path:
        filename, extension = os.path.splitext(path)
        if not extension or extension != '.nii':
            extension = '.nii'
            path = filename + extension
            nib.nifti1.save(nii, path)
    return nii

# %%
