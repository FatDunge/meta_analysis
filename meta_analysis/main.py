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
import copy

import numpy as np
import pandas as pd

from . import meta
from . import data
from . import utils
from . import mask

def load_centers_data(center_dict):
    """ generate list of Center instance
    Args:
        center_dict: dict of dict of filepathes.
                     {center1_name:{group1_label:[filepath1, filepath2],
                                   group2_label:[filepath3, filepath4],
                                   ...}
                      center2_name:{...}}
    Return:
        center_dict: dict of dict of datas.
                     {center1_name:{group1_label:[data1, data2],
                                   group2_label:[data3, data4],
                                   ...}
                      center2_name:{...}}
    """
    for center_name, group_dict in center_dict.items():
        for label, filepathes in group_dict.items():
            array = utils.load_arrays(filepathes)
            center_dict[center_name][label] = array
    return center_dict

def voxelwise_meta_analysis(center_dict, label1, label2,
                            _mask=None, is_filepath=True,
                            model='random', method='cohen_d'):
    """ perform voxelwise meta analysis
    Args:
        center_dict: dict of dict of group filepathes. pass to gen_array_center()
                    {center1:{group1:[filepath1, filepath2],
                              group2:[filepath3, filepath4]}
                     center2:{...}}
        label1: label of experimental group
        label2: label of control group
        _mask: Mask instance, use to mask array, will only caculate mask region.
        is_filepath: bool, is filepath or ndarray. pass to gen_array_center()
        model: 'fixed' or 'random', meta analysis model.
        method: str, ways to caculate effect size
    Return:
        results: ndarray, return from VoxelMeta.main()
    """
    if is_filepath:
        center_dict = load_centers_data(center_dict)
    
    origin_shape = None
    flatten_shape = None
    for k, group_dict in center_dict.items():
        for label, datas in group_dict.items():
            if origin_shape is None:
                origin_shape = datas[0].shape
            datas = np.reshape(datas, (len(datas), -1))
            center_dict[k][label] = datas
            if flatten_shape is None:
                flatten_shape = datas[0].shape
    if _mask is not None:
        if _mask.get_shape() == origin_shape:
            _mask = mask.Mask(_mask.data.flatten())
        elif _mask.get_shape() == flatten_shape:
            pass
        else:
            raise AssertionError('Mask shape couldn\'t fit with data')
    
    center_mean_dict = copy.deepcopy(center_dict)
    center_std_dict = copy.deepcopy(center_dict)
    center_count_dict = copy.deepcopy(center_dict)
    for center_name, group_dict in center_dict.items():
        for label, datas in group_dict.items():
            mean, std, count = utils.cal_mean_std_n(datas)
            center_mean_dict[center_name][label] = mean
            center_std_dict[center_name][label] = std
            center_count_dict[center_name][label] = count

    if _mask is not None:
        indexes = _mask.get_nonzero_index()
    else:
        indexes = np.transpose(np.nonzero(np.ones(flatten_shape)))

    m = meta.Meta()
    results_array = None
    import time
    st = time.time()
    for index in indexes:
        center_list = []
        for center_name, group_dict in center_dict.items():
            groups = []
            for label, _ in group_dict.items():
                mean = center_mean_dict[center_name][label][index][0]
                std = center_std_dict[center_name][label][index][0]
                count = center_count_dict[center_name][label]
                group = data.Group(label=label, mean=mean, std=std, count=count)
                groups.append(group)
            center = data.Center(center_name, groups)
            center_list.append(center)
        centers = data.Centers(center_list)
        studies = centers.gen_studies(label1, label2, method)
        result_model = m.caculate(studies, model)
        results = result_model.get_results()
        if results_array is None:
            results_len = len(results)
            results_array = np.zeros(flatten_shape+(results_len,))
        results_array[index] = results
    results_array = np.transpose(results_array)
    results_array = np.reshape(results_array, (results_len,)+origin_shape)
    et = time.time()
    els = et - st
    return results_array

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
    """
    pass

def csv_meta_analysis(csvpath, header=0, method='cohen_d', model='random'):
    """ perform meta analysis based on csv file
    Args:
        csvpath: csv filepath,
                 csv example:
                    center_name, m1, s1, n1, m2, s2, n2
                    center1, 1, 1, 10, 2, 2, 20
                    center2, 1.2, 2, 15, 2.2, 2, 15
        header: 0 or None, pandas.read_csv args.
                None means no header
        model: 'fixed' or 'random', meta analysis model.
        method: str, ways to caculate effect size
    Return:
        results: Model instance
    """
    df = pd.read_csv(csvpath, header=header, index_col=0)
    studies = []
    for index, row in df.iterrows():
        m1,s1,n1,m2,s2,n2 = row.values
        group1 = data.Group(0, mean=m1, std=s1, count=n1)
        group2 = data.Group(1, mean=m2, std=s2, count=n2)
        center = data.Center(index, [group1, group2])
        studies.append(center.gen_study(0,1,method))
    _meta = meta.Meta()
    result_model = _meta.caculate(studies, _model=model)
    return result_model

