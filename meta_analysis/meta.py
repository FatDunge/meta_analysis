import nibabel as nib
import numpy as np

import meta_analysis.model

class Meta(object):
    def __init__(self):
        super().__init__()
    
    def load_data(self):
        pass

    def caculate(self, studies, model, index=None):
        total_weight = 0
        effect_sizes = []
        variances = []
        for study in studies:
            es = study.get_effect_size(index)
            v = study.get_variance(index)

            effect_sizes.append(es)
            variances.append(v)

        if model == 'Fixed':
            m = meta_analysis.model.FixedModel(effect_sizes, variances)
        elif model == 'Random':
            m = meta_analysis.model.RandomModel(effect_sizes, variances)
        return m

class VoxelMeta(Meta):
    def __init__(self, centers):
        super().__init__()
        self.centers = centers
        self.check_centers()
        self.results_array = None

    def check_centers(self):
        shape = self.centers[0].shape
        for center in self.centers:
            assert shape == center.shape

    def main(self, label1, label2, model, method, mask=None):
        array_shape = self.centers[0].shape
        if not mask:
            mask = np.ones(array_shape)
        assert mask.shape == array_shape
        i = 0
        for index, x in np.ndenumerate(mask):
            if x > 0:
                studies = []
                for center in self.centers:
                    study = center.gen_study(label1, label2, method)
                    studies.append(study)
                m = self.caculate(studies, model, index)
                results = m.get_results()
                if i == 0:
                    results_array = np.zeros(((len(results),) + array_shape))
                for j in range(len(results)):
                    results_array[j][index] = results[j]
        self.results_array = results_array

