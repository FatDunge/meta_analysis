import numpy as np
import scipy

class Model(object):
    def __init__(self):
        super().__init__()

    def caculate(self):
        raise NotImplementedError

def inverse_variance(variance):
    return 1 / variance

def get_confidence_intervals(effect_size, standard_error):
    lower_limits = effect_size - 1.96 * standard_error
    upper_limits = effect_size + 1.96 * standard_error
    return lower_limits, upper_limits

def get_heterogeneity(effect_sizes, total_es, weights):
    diff_es_square = (effect_sizes - total_es) ** 2
    q = np.sum(np.multiply(weights, diff_es_square))
    return q

def get_z_value(total_es, standard_error, x0=0):
    return (total_es - x0) / standard_error

def get_p_from_z(z, one_side=False):
    if one_side:
        p_values = scipy.stats.norm.sf(abs(z))
    else:
        p_values = scipy.stats.norm.sf(abs(z))*2

class FixedModel(Model):
    def __init__(self):
        super().__init__()

    def caculate(self, effect_sizes, variances):
        weights = np.reciprocal(variances)
        total_es = np.sum(np.multiply(effect_sizes, weights)) /\
                   np.sum(weights)
        total_v = 1 / np.sum(weights)
        standard_error = np.sqrt(total_v)

        ll, ul = get_confidence_intervals(total_es, standard_error)
        q = get_heterogeneity(effect_sizes, total_es, weights)
        z = get_z_value(total_es, standard_error)
        p = get_p_from_z(z)



    
