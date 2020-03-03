import numpy as np
from scipy.stats import norm

def inverse_variance(variance):
    return 1 / variance

def get_confidence_intervals(effect_size, standard_error):
    lower_limits = effect_size - 1.96 * standard_error
    upper_limits = effect_size + 1.96 * standard_error
    return lower_limits, upper_limits

def get_heterogeneity(effect_sizes, total_effect_size, weights):
    diff_es_square = (effect_sizes - total_effect_size) ** 2
    q = np.sum(np.multiply(weights, diff_es_square))
    return q

def get_z_value(total_effect_size, standard_error, x0=0):
    return (total_effect_size - x0) / standard_error

def get_p_from_z(z, one_side=False):
    if one_side:
        p_value = norm.sf(abs(z))
    else:
        p_value = norm.sf(abs(z)) * 2
    return p_value

class Model(object):
    def __init__(self, effect_sizes, variances):
        super().__init__()
        self.effect_sizes = np.asarray(effect_sizes)
        self.variances = np.asarray(variances)
        self.gen_weights()

        self.caculate()
    
    def gen_weights(self):
        raise NotImplementedError("Generate Weight Method Not Implemented")

    def caculate(self):
        effect_sizes = self.effect_sizes
        variances = self.variances
        weights = self.weights

        total_effect_size = np.sum(np.multiply(effect_sizes, weights)) /\
                   np.sum(weights)
        total_variance = 1 / np.sum(weights)
        standard_error = np.sqrt(total_variance)

        lower_limits, upper_limits = get_confidence_intervals(total_effect_size, standard_error)
        q = get_heterogeneity(effect_sizes, total_effect_size, weights)
        z = get_z_value(total_effect_size, standard_error)
        p = get_p_from_z(z)

        self.total_effect_size = total_effect_size
        self.total_variance = total_variance
        self.standard_error = standard_error
        self.lower_limits = lower_limits
        self.upper_limits = upper_limits
        self.q = q
        self.z = z
        self.p = p

    def get_results(self):
        return (self.total_effect_size,
                self.total_variance, self.standard_error,
                self.lower_limits, self.upper_limits,
                self.q, self.z, self.p)

class FixedModel(Model):
    def __init__(self, effect_sizes, variances):
        super().__init__(effect_sizes, variances)

    def gen_weights(self):
        effect_sizes = self.effect_sizes
        variances = self.variances

        self.weights = np.reciprocal(variances)

class RandomModel(Model):
    def __init__(self, effect_sizes, variances):
        super().__init__(effect_sizes, variances)

    def gen_weights(self):
        effect_sizes = self.effect_sizes
        variances = self.variances

        fixed_weights = np.reciprocal(variances)
        mean_effect_size = np.mean(effect_sizes)
        Q = np.sum(np.square(variances-mean_effect_size)/variances)
        df = len(variances) - 1
        C = np.sum(fixed_weights) - np.sum(np.square(fixed_weights)) / np.sum(fixed_weights)
        tau_square = (Q - df) / C
        if tau_square < 0:
            tau_square = 0

        self.weights = np.reciprocal(variances + tau_square)