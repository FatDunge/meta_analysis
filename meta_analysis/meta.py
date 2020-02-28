import nibabel as nib
import numpy as np

import meta_analysis.model as model

class Meta(object):
    def __init__(self):
        super().__init__()
    
    def load_data(self):
        pass

    def caculate(self, studies, model, index):
        total_weight = 0
        effect_sizes = []
        variances = []
        for study in studies:
            es = study.get_effect_size(index)
            v = study.get_variance()

            effect_sizes.append(es)
            variances.append(v)

        m = model.FixedModel(effect_sizes, variances)
        results = m.get_results()

class VoxelMeta(Meta):
    def __init__(self, centers):
        super().__init__()
        self.centers = centers

    def main(self, label1, label2, model, method, mask=None):
        if not mask:
            mask = np.zeros(centers[0].array_shape)
        for index, x in np.ndenumerate(mask):
            studies = []
            for center in centers:
                study = center.gen_study(label1, label2, method)
                studies.append(study)
            super.caculate(studies, index)
