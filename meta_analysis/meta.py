import nibabel as nib
import numpy as np

class Meta(object):
    def __init__(self, model, effect_size):
        super().__init__()
        self.model = model
        self.effect_size = effect_size
    
    def load_data(self):
        pass

    def caculate(self, studies):
        total_weight = 0
        for study in studies:
            if self.effect_size == 'cohen_d':
                es = study.cohen_d()

def get_mean_std_n(arrays):
    mean = np.mean(arrays, axis=0)
    std = np.std(arrays, axis=0)
    n = arrays.shape[0]
    return mean, std, n

class VoxelMeta(Meta):
    def __init__(self, centers):
        super().__init__()
        self.centers = centers

    def main(self, mask=None):
        if not mask:
            tmp = np.zeros(centers[0].array_shape)
            for index, x in np.ndenumerate(tmp):
                studies = []
                for center in centers:
                    study = center.gen_study(index)
                    studies.append(study)
                super.caculate(studies)
