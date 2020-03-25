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
import matplotlib.pyplot as plt
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
        lower_limits = []
        upper_limits = []
        for study in studies:
            es = study.get_effect_size()
            v = study.get_variance()
            ll, ul = study.get_confidence_intervals()
            effect_sizes.append(es)
            variances.append(v)
            lower_limits.append(ll)
            upper_limits.append(ul)

        self.effect_sizes = np.asarray(effect_sizes)
        self.variances = np.asarray(variances)
        self.lower_limits = np.asarray(lower_limits)
        self.upper_limits = np.asarray(upper_limits)
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
        total_standard_error = np.sqrt(total_variance)

        total_lower_limit, total_upper_limit = get_confidence_intervals(total_effect_size, total_standard_error)
        q = get_heterogeneity(effect_sizes, total_effect_size, weights)
        z = get_z_value(total_effect_size, total_standard_error)
        p = get_p_from_z(z)

        self.total_effect_size = total_effect_size
        self.total_variance = total_variance
        self.total_standard_error = total_standard_error
        self.total_lower_limit = total_lower_limit
        self.total_upper_limit = total_upper_limit
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
        dpi = 100
        width = 14
        height = 2+len(self.studies)
        grid_width = 1 / width
        grid_height = 1 / height
        fig, ax = plt.subplots(figsize=(width, height),
                               dpi=dpi)
        # set default to invisible
        ax.tick_params(left=False, bottom=False, labelbottom=False, labelleft=False)
        for spine in ax.spines.values():
            spine.set_visible(False)
        
        x, y = 0, 0
        ha = 'center'
        #draw header
        study_x = x
        experimental_group_x = x + grid_width * 3
        control_group_x = x + grid_width * 6
        froest_plot_x = x + grid_width * 9
        effect_size_x = x + grid_width * 11
        lower_limit_x = x + grid_width * 12
        upper_limit_x = x + grid_width * 13
        interval_x = (lower_limit_x + upper_limit_x) / 2
        weight_x = x + grid_width * 14
        header_y = 1

        ax.text(study_x, 1, 'Study')
        ax.text(experimental_group_x, header_y, 'experimental group',ha=ha)
        ax.text(control_group_x, header_y, 'control group',ha=ha)
        ax.text(froest_plot_x, header_y, 'forest plot',ha=ha)
        ax.text(effect_size_x, header_y, 'effect size',ha=ha)
        ax.text(interval_x, header_y, '95% interval',ha=ha)
        ax.text(weight_x, header_y, 'weight',ha=ha)
        #draw subheader
        subheader_y = 1 - grid_height / 2
        eg_mean_x = experimental_group_x - grid_width
        eg_std_x = experimental_group_x
        eg_count_x = experimental_group_x + grid_width
        cg_mean_x = control_group_x - grid_width
        cg_std_x = control_group_x
        cg_count_x = control_group_x + grid_width

        ax.text(eg_mean_x, subheader_y, 'mean',ha=ha)
        ax.text(eg_std_x, subheader_y, 'std',ha=ha)
        ax.text(eg_count_x, subheader_y, 'count',ha=ha)
        ax.text(cg_mean_x, subheader_y, 'mean',ha=ha)
        ax.text(cg_std_x, subheader_y, 'std',ha=ha)
        ax.text(cg_count_x, subheader_y, 'count',ha=ha)

        total_eg_count = 0
        total_cg_count = 0

        row_y = 1 - grid_height
        ax.axhline((subheader_y+row_y)/2, color='black')
        for i, (study, effect_size, weight,
                lower_limit, upper_limit) in enumerate(
                    zip(self.studies, self.effect_sizes, self.weights,
                        self.lower_limits, self.upper_limits), 1):
            row_y = 1-grid_height*i
            ax.text(study_x, row_y, study.name)
            eg_mean, eg_std, eg_count = study.group_experimental.get_mean_std_count()
            cg_mean, cg_std, cg_count = study.group_control.get_mean_std_count()
            ax.text(eg_mean_x, row_y, eg_mean)
            ax.text(eg_std_x, row_y, eg_std)
            ax.text(eg_count_x, row_y, eg_count)
            ax.text(cg_mean_x, row_y, cg_mean)
            ax.text(cg_std_x, row_y, cg_std)
            ax.text(cg_count_x, row_y, cg_count)
            ax.text(effect_size_x, row_y, '{:.2f}'.format(effect_size), ha=ha)
            ax.text(lower_limit_x, row_y, '[{:.2f}'.format(lower_limit), ha=ha)
            ax.text(upper_limit_x, row_y, '{:.2f}]'.format(upper_limit), ha=ha)
            ax.text(weight_x, row_y, '{:.2f}%'.format(weight*100), ha=ha)

            total_eg_count += eg_count
            total_cg_count += cg_count
        # draw total
        total_y = grid_height
        ax.axhline((row_y+total_y)/2, color='black')
        ax.text(study_x, total_y, 'Total')
        ax.text(eg_count_x, total_y, total_eg_count)
        ax.text(cg_count_x, total_y, total_cg_count)
        ax.text(effect_size_x, total_y, '{:.2f}'.format(self.total_effect_size), ha=ha)
        ax.text(lower_limit_x, total_y, '[{:.2f}'.format(self.total_lower_limit), ha=ha)
        ax.text(upper_limit_x, total_y, '{:.2f}]'.format(self.total_upper_limit), ha=ha)
        ax.text(weight_x, total_y, '{:.2f}%'.format(np.sum(self.weights)*100), ha=ha)
        
        # draw summary
        summary_y = grid_height / 2
        ax.text(0, summary_y, 'Heterogeneity:{:.2f}'.format(self.q), ha=ha)
        ax.text(2*grid_width, summary_y, 'z-value:{:.2f}'.format(self.z), ha=ha)
        ax.text(4*grid_width, summary_y, 'p-value:{:.2f}'.format(self.p), ha=ha)
        plt.show()


class FixedModel(Model):
    def __init__(self, studies):
        super().__init__(studies)

    def gen_weights(self):
        weights = np.reciprocal(self.variances)
        self.weights = weights / np.sum(weights)

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

        weights = np.reciprocal(variances + tau_square)
        self.weights = weights / np.sum(weights)