import nibabel as nib
import numpy as np

import meta_analysis.model as model

class Meta(object):
    def __init__(self):
        super().__init__()
    
    def load_data(self):
        pass

    def caculate(self, studies, model):
        total_weight = 0
        effect_sizes = []
        variances = []
        for study in studies:
            es = study.get_effect_size()
            v = study.get_variance()

            effect_sizes.append(es)
            variances.append(v)

        m = model.FixedModel(effect_sizes, variances)
        results = m.get_results()


def get_mean_std_n(arrays):
    mean = np.mean(arrays, axis=0)
    std = np.std(arrays, axis=0)
    n = arrays.shape[0]
    return mean, std, n

class VoxelMeta(Meta):
    def __init__(self, centers):
        super().__init__()
        self.centers = centers

    def main(self, model, method, mask=None):
        if not mask:
            tmp = np.zeros(centers[0].array_shape)
            for index, x in np.ndenumerate(tmp):
                studies = []
                for center in centers:
                    study = center.gen_study(index, method)
                    studies.append(study)
                super.caculate(studies)
