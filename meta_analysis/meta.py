import nibabel as nib
import numpy as np

import meta_analysis.model

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

        if model.lower() == 'fixed':
            m = meta_analysis.model.FixedModel(effect_sizes, variances)
        elif model.lower() == 'random':
            m = meta_analysis.model.RandomModel(effect_sizes, variances)
        return m

class VoxelMeta(Meta):
    def __init__(self, centers):
        super().__init__()
        self.centers = centers
        self.results_array = None

    def check_centers(self):
        base_center = self.centers[0]
        for center in self.centers:
            assert base_center.is_same_shape(center)

    def main(self, label1, label2, model, method, mask=None):
        array_shape = self.centers[0].shape
        if mask is None:
            mask = np.ones(array_shape)
        assert mask.shape == array_shape
        i = 0
        for index, x in np.ndenumerate(mask):
            if x > 0:
                studies = []
                for center in self.centers:
                    study = center.gen_study(label1, label2, method, index)
                    studies.append(study)
                m = super().caculate(studies, model)
                results = m.get_results()
                if i == 0:
                    results_array = np.zeros(((len(results),) + array_shape))
                for j in range(len(results)):
                    results_array[j][index] = results[j]
        self.results_array = results_array
        return self.results_array

