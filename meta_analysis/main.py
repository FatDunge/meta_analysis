from . import meta
from . import data

def gen_array_center(center_dict, is_filepath):
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
    centers = gen_array_center(center_dict, is_filepath)
    
    m = meta.VoxelMeta(centers)
    results = m.main(label1, label2, model, method, mask)
    return results

def region_volume_meta_analysis(center_dict, label1, label2, 
                                mask, is_filepath=True,
                                model='random', method='cohen_d'):
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
    return result
