""" main file of full process, from load data to get results

Function:
    gen_array_center(center_dict, is_filepath): generate center contains ArrayGroup from dict
    voxelwise_meta_analysis(center_dict, label1, label2,
                            mask, is_filepath, model, method): perform voxelwise meta analysis
    region_volume_meta_analysis(center_dict, label1, label2, 
                            mask, is_filepath, model, method): perform region volume meta analysis

Author: Kang Xiaopeng
Data: 2020/03/05
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

from . import meta
from . import data

def gen_array_center(center_dict, is_filepath):
    """ generate list of Center instance
    Args:
        center_dict: dict of dict of group filepathes.
                    {center1:{group1:[filepath1, filepath2],
                              group2:[filepath3, filepath4]}
                     center2:{...}}
        is_filepath: bool, is filepath or ndarray
    Return:
        centers: list of Center instance
    """
    centers = []
    for center_name, group_dict in center_dict.items():
        groups = []
        for label, datas in group_dict.items():
            if is_filepath:
                group = data.ArrayGroup(label, datapathes=datas)
            else:
                group = data.ArrayGroup(label, datas=datas)
            groups.append(group)
        center = data.Center(center_name, groups)
        centers.append(center)
    return centers

def voxelwise_meta_analysis(center_dict, label1, label2,
                            mask=None, is_filepath=True,
                            model='random', method='cohen_d'):
    """ perform voxelwise meta analysis
    Args:
        center_dict: dict of dict of group filepathes. pass to gen_array_center()
                    {center1:{group1:[filepath1, filepath2],
                              group2:[filepath3, filepath4]}
                     center2:{...}}
        label1: label of experimental group
        label2: label of control group
        mask: Mask instance, use to mask array, will only caculate mask region.
        is_filepath: bool, is filepath or ndarray. pass to gen_array_center()
        model: 'fixed' or 'random', meta analysis model.
        method: str, ways to caculate effect size
    Return:
        results: ndarray, return from Meta.main()
    """
    
    centers = gen_array_center(center_dict, is_filepath)
    
    m = meta.VoxelMeta(centers)
    results = m.main(label1, label2, model, method, mask)
    return results

def region_volume_meta_analysis(center_dict, label1, label2, 
                                mask, is_filepath=True,
                                model='random', method='cohen_d'):
    """ perform region volume meta analysis
    Args:
        center_dict: dict of dict of group filepathes. pass to gen_array_center()
                    {center1:{group1:[filepath1, filepath2],
                              group2:[filepath3, filepath4]}
                     center2:{...}}
        label1: label of experimental group
        label2: label of control group
        mask: Mask instance, use to mask array, will only caculate mask region.
        is_filepath: bool, is filepath or ndarray. pass to gen_array_center()
        model: 'fixed' or 'random', meta analysis model.
        method: str, ways to caculate effect size
    Return:
        results: dict of tuple, {region_label1: result1, ...}
    """
    centers = gen_array_center(center_dict, is_filepath)

    region_labels = mask.get_labels()
    results = {}
    for region_label in region_labels:
        studies = []
        for center in centers:
            study = center.gen_region_study(label1, label2, method, mask, region_label)
            studies.append(study)
        m = meta.Meta()
        result = m.caculate(studies, model)
        results[region_label] = result
    return results
