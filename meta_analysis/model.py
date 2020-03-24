""" model module as meta analysis model

Function:
    inverse_variance(variance): inverse variance
    get_confidence_intervals(effect_size, standard_error): caculate 95% confidence intervals
    get_heterogeneity(effect_sizes, total_effect_size, weights)： caculate heterogeneity
    get_z_value(total_effect_size, standard_error, x0): caculate z test value
    get_p_from_z(z, one_side): caculate p value from z value

Class:
    Model(object): basic model
    FixedModel(Model): Fixed model, assume one true effect size that underlies all the studies.
    RandomModel(Model): Random model, assume effect size varies from study to study.

Author: Kang Xiaopeng
Data: 2020/02/21
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

import numpy as np
from scipy.stats import norm

def inverse_variance(variance):
    return 1 / variance

def get_confidence_intervals(effect_size, standard_error):
    # 95% confidence intervals
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
    """Basic class to perform check effect size from all study.

    Attributes:
        effect_sizes: list, all studies' effect size
        variances: list, all studies' variance
        weights: list, all studies' weight
        total_effect_size: float, combined effect size
        total_variance: float, combined variance 
        standard_error: float, standard error of meta analysis
        lower_limits: float, 95% lower confidence intervals
        upper_limits: float, 95% upper confidence intervals
        q: float, heterogeneity
        z: float, z test value
        p: float, p value

    Function:
        gen_weights(): caculate weight from effect_sizes and variances
        caculate(): caculate results
    """

    def __init__(self, studies):
        super().__init__()
        self.studies = studies
        effect_sizes = []
        variances = []

        for study in studies:
            es = study.get_effect_size()
            v = study.get_variance()

            effect_sizes.append(es)
            variances.append(v)

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

    def plot_forest(self, show_weights=True, show_effect_sizes=True,
                    sort=False):
        pass


class FixedModel(Model):
    def __init__(self, studies):
        super().__init__(studies)

    def gen_weights(self):
        self.weights = np.reciprocal(self.variances)

class RandomModel(Model):
    def __init__(self, studies):
        super().__init__(studies)

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