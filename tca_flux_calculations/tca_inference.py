"""

"""

import numpy as np
import pandas as pd

import matplotlib

font = {'family': 'normal',
        'weight': 'normal',
        'size': 16}

matplotlib.rc('font', **font)
import matplotlib.pyplot as plt

from scipy.optimize import least_squares
from scikits.odes.odeint import odeint
from scipy.stats import t as t_distribtion

"""
Parameter sections with helpful deifnitions
"""

DEFAULT_MAX_FLUX = 1e5

extra_options = {'old_api': False}

BARTMAN_DATA_MAPPING = {'lactate': 'lac', 'aspartate': 'asp', 'malate': 'mal',
                        'succinate': 'succ', 'glutamate': 'glu',
                        'citrate/isocitrate': 'cit', 'citrate': 'cit',
                        'a-ketoglutarate': 'akg', 'acetylcoa': 'accoa', 'oxolacetate': 'oaa',
                        'pyruvate': 'pyr'}


class InstatFluxFitter(object):
    """
    This class defines a set of methods required to infer fluxes from instationary labeling data
    """

    def __init__(self,
                 tca_model,
                 metabolite_ix_hash_map,
                 parameters,
                 data_model_mapping=BARTMAN_DATA_MAPPING,
                 max_flux=DEFAULT_MAX_FLUX,
                 M=10, ):

        self.tca_model = tca_model

        self.metablite_ix_hash_map = metabolite_ix_hash_map

        self.data_model_mapping = data_model_mapping

        self.max_flux = max_flux

        self.parameters = parameters
        self.parameter_mapping = {k: i for i, k in enumerate(parameters)}

        self.no_lablel = [k for k in metabolite_ix_hash_map if k.endswith('_0b0')]

        isotopomer_map = {k: (k.split('_')[0], sum([int(x) for x in k.split('_')[1][2:]]))
                          for k in metabolite_ix_hash_map}

        self.isotopoimer_hash_map = {k: {i: [] for i in range(0, M)} for k in data_model_mapping.values()}
        for isomer in metabolite_ix_hash_map:
            ix = metabolite_ix_hash_map[isomer]
            species, label = isotopomer_map[isomer]
            self.isotopoimer_hash_map[species][label].append(ix)

    def pool_sizes_from_data(self, data, small_pool=0.1):
        """

        :param small_pool:
        :return:
        """

        pool_sizes = dict()

        for compound_name, row in data.iterrows():
            comp = self.data_model_mapping[compound_name]
            pool_size = -1
            while pool_size < 0:
                pool_size = np.random.randn() * row.sdconc + row.meanconc

            pool_sizes[comp] = pool_size

        min_poolsize = min(pool_sizes.values())

        for met in self.data_model_mapping.values():
            if not met in pool_sizes:
                pool_sizes[met] = min_poolsize * small_pool

        return pool_sizes

    def get_pools_and_params(self, fitting_result):
        """

        :param fitting_result: pd.Series
        :return:
        """
        parameters = np.array([fitting_result[p] for p in self.parameters])
        pool_sizes = {p: fitting_result[p] for p in self.data_model_mapping.values()}
        return parameters, pool_sizes

    def solve_odes(self, time, params, pool_sizes):
        """

        :param params:
        :param pool_sizes:
        :return:
        """
        L0 = np.zeros(len(self.metablite_ix_hash_map))
        L0[np.array([self.metablite_ix_hash_map[i] for i in self.no_lablel])] = 1.0

        func = lambda t, y, dydt: self.tca_model(t, y, dydt, params, pool_sizes)

        sol = odeint(func, time, L0,
                     max_steps=1e9,
                     atol=1e-6,
                     rtol=1e-6,
                     )
        return sol

    def residuals(self, parameters, data, pool_sizes, ):

        timepoints = sorted(data.time.unique())
        timepoints = np.append(0, timepoints)
        timepoints = np.append(timepoints, 200)
        sol = self.solve_odes(timepoints, parameters, pool_sizes, )

        time = sol.values.t
        values = sol.values.y

        if len(time) != len(timepoints):
            print('Solver error ')
            return 1e3

        res = []

        time_values = pd.DataFrame(values, index=time)
        unique = data.groupby(['Compound', 'C_Label']).size().reset_index().rename(columns={0: 'count'})

        for i, row in unique.iterrows():
            try:
                if row['Compound'] == 'lactate':
                    raise KeyError('No lactate pls!')

                indices_cumomers = self.isotopoimer_hash_map[self.data_model_mapping[row['Compound']]][row['C_Label']]

                data_slice = data[(data['Compound'] == row['Compound'])
                                  & (data['C_Label'] == row['C_Label'])]

                cumomer_labeling = time_values.loc[data_slice['time'], indices_cumomers].sum(axis=1)

                r = (data_slice['fraction'].values - cumomer_labeling.values)

                res.append(r)


            except KeyError:
                pass

        res_lsq = np.concatenate(res)

        return res_lsq

    def plot_data_fit(self, sol, data, folder, name, M=3, ):

        compounds = data.Compound.unique()
        compounds = [c for c in compounds if c != 'lactate']

        fig, axes = plt.subplots(nrows=3, ncols=len(compounds), sharex=True, sharey=True, figsize=(10, 7))

        for i, compound in enumerate(compounds):
            for j, label in enumerate(range(1, M + 1)):
                this_axis = axes[j, i]

                indices_cumomers = self.isotopoimer_hash_map[self.data_model_mapping[compound]][label]

                # steady_state_labeling_metabolite = sol.values.y[-1,indices_metabolite].sum()
                steady_state_labeling_metabolite = 1.0

                data_slice = data[(data['Compound'] == compound)
                                  & (data['C_Label'] == label)]
                if not type(sol) == list:

                    this_axis.plot(sol.values.t, sol.values.y[:, indices_cumomers].sum(axis=1)
                                   / steady_state_labeling_metabolite, 'k--')
                else:
                    for s in sol:
                        this_axis.plot(s.values.t, s.values.y[:, indices_cumomers].sum(axis=1)
                                       / steady_state_labeling_metabolite, 'k-', alpha=0.5)

                # Plot data
                this_axis.scatter(data_slice.time, data_slice.fraction / steady_state_labeling_metabolite, c='k')

        pad = 5  # in points

        for ax, col in zip(axes[0, :], compounds):
            ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
                        xycoords='axes fraction', textcoords='offset points',
                        size='large', ha='center', va='baseline')

        plt.tight_layout()
        plt.savefig(folder + '/' + name + '.png')
        plt.close()

    def export_ode_solution(self, sol, folder='output', name='tissue', ):
        compounds = self.data_model_mapping.keys()
        data = {'time': sol.values.t}
        for i, compound in enumerate(compounds):
            for j, label in enumerate(range(1, 4)):
                indices_cumomers = self.isotopoimer_hash_map[self.data_model_mapping[compound]][label]
                data[compound + ' M+' + str(label)] = sol.values.y[:, indices_cumomers].sum(axis=1)

        df = pd.DataFrame(data=data)
        df.to_csv(folder + '/' + name + '.csv')

    def fit_data(self, fn_inpt, verbose=0):

        data, pool_sizes, P0, param_bounds = fn_inpt

        lsq_parameter_bounds = list(zip(*param_bounds))
        res = least_squares(self.residuals, P0,
                            args=(data, pool_sizes,),
                            bounds=lsq_parameter_bounds,
                            verbose=verbose,
                            )

        return res, pool_sizes

    def make_input_data(self, tissue, pool_size_data, labeling_data, tissue_time_scale, N_init, 
                        M=3, value_column="conc_micromolar"):

        # Select data for one tissue
        tissue_for_pool = tissue
        ps_data = pool_size_data[pool_size_data['tissue'] == tissue_for_pool]

        ps_data = pd.concat([ps_data.groupby('Compound')[value_column].mean(),
                             ps_data.groupby('Compound')[value_column].std()], axis=1)

        ps_data.columns = ['meanconc', 'sdconc']

        selected_data = labeling_data[(labeling_data['tissue'] == tissue)
                                      & (labeling_data['omitFromFit'] == False)
                                      & (labeling_data['C_Label'] <= M)]

        upper_flux_guess = ps_data[ps_data.index != 'lactate'].sum()['meanconc'] \
                           * tissue_time_scale[tissue] * 2.0

        param_bounds = PARAMETER_BOUNDS(self.max_flux, 1000)
        init_param_bounds = PARAMETER_BOUNDS(upper_flux_guess, 10, )

        P0 = initial_parameters(N_init, init_param_bounds)

        pool_sizes_input = []
        for i in range(N_init):
            this_pool_sizes = self.pool_sizes_from_data(ps_data)
            pool_sizes_input.append(this_pool_sizes)

        # data,pool_sizes,P0,param_bounds,VMAX):
        function_inputs = [(selected_data, pools, p0, param_bounds)
                           for p0, pools in zip(P0, pool_sizes_input)]

        return function_inputs

    def export_result(self, results, tissue=None):
        """
        Export raw output to table format

        :param results: A list of fitting results (outputs of self.fit_data)
        :param tissue: Tissue name for annotation
        :return: pd.DataFrame [ Pool concentrations ,Parameters,  cost function, tissue (optional)]
        """

        this_results = []
        for res in results:
            # Cacl R2
            pool_sizes = res[1]
            parameters = {r: res[0]['x'][v] for r, v in self.parameter_mapping.items()}
            cost = {'cost': res[0]['cost']}
            row_dict = dict()
            row_dict.update(pool_sizes)
            row_dict.update(parameters)
            row_dict.update(cost)
            if not tissue is None:
                row_dict.update({'tissue': tissue})

            row = pd.DataFrame(row_dict, index=[0, ])
            this_results.append(row)

        this_results = pd.concat(this_results).reset_index(drop=True)

        return this_results

    # Boot straping
    def boostrap_parameters(self,
                            result_dataset,
                            k=20, n=100,
                            columns=['TCA'],
                            alpha=5.0, ):

        # Probably filter data set based on top X
        this_dataset = result_dataset.sort_values('cost')[:k]
        this_dataset = this_dataset.reset_index(drop=True)

        if this_dataset.shape[0] < 2:
            raise RuntimeError()

        scores = list()

        m = this_dataset.shape[0]
        for _ in range(n):
            # bootstrap sample
            indices = np.random.randint(0, m, m)
            sample = this_dataset.loc[indices]

            # calculate and store statistic
            statistic = np.median(sample[columns], axis=0)

            scores.append(statistic)

        lower_p = alpha / 2.0
        upper_p = (100 - alpha) + (alpha / 2.0)

        lower = np.percentile(scores, lower_p, axis=0)
        upper = np.percentile(scores, upper_p, axis=0)
        median = np.percentile(scores, 50, axis=0)

        return lower, upper, median,

    def percentile_boostrap_p_value(self,
                                    result1, result2,
                                    k=20, n=1000,
                                    columns=['TCA'], ):

        this_dataset_1 = result1.sort_values('cost')[:k].reset_index()
        this_dataset_2 = result2.sort_values('cost')[:k].reset_index()

        diff = this_dataset_1[columns] - this_dataset_2[columns]
        t = diff.median()
        t_star = list()

        m = diff.shape[0]
        for _ in range(n):
            # bootstrap sample
            indices = np.random.randint(0, m, m)
            sample = diff.loc[indices] - t

            # calculate and store statistic
            statistic = np.median(sample[columns], axis=0)

            t_star.append(statistic)

        t_stat = abs(np.median(t_star) - t) / np.std(t_star)

        # Two-sided t-testnp.median(t_star)
        p_value = (1 - t_distribtion.cdf(x=t_stat, df=n - 1)) * 2.0

        return p_value

    def find_confidence_intervals(self, results, columns=['TCA', ], alpha=5.0, k=20, n=100, ):
        """
        Determine the confidence intervall of the output parameters

        :param results: List of results from sel.fit_data()
        :param columns: List of parameters
        :param alpha: Confidence level
        :param k: Top k fits to bootstrap
        :param n: Number of resamples for the bootstrapping
        :return: dict("parameter" : dict() ) median and ci bounds for each parameter defined in columns
        """

        lower, upper, median, = self.boostrap_parameters(results,
                                                         columns=columns,
                                                         alpha=alpha,
                                                         k=k,
                                                         n=n, )

        confidence_intervals = dict()

        for i, c in enumerate(columns):
            confidence_intervals[c] = {'median': median[i], 'ci_upper': upper[i], 'ci_lower': lower[i]}

        return confidence_intervals


"""
Ultility functions
"""


def PARAMETER_BOUNDS(VMAX, RMAX, zero=0, r=(0, 1), da=(0, 1), dg=(-1, 1)):
    return [dg, (zero, VMAX), (zero, VMAX), (zero, VMAX),  #
            (1, RMAX), (1, RMAX), (1, RMAX), (1, RMAX),
            (1, RMAX), (1, RMAX), (1, RMAX), (zero, 1), da, r]


def initial_parameters(M, param_bounds):
    # RANDOM = lhsmdu.sample(len(param_bounds),M).T
    RANDOM = np.random.rand(len(param_bounds), M).T
    initial_params = []
    for R in np.asarray(RANDOM):
        p0 = [r * (u - l) + l for r, (l, u) in zip(R, param_bounds)]
        initial_params.append(p0)
    return initial_params


def ci_table(tissue_cis, parameter):
    table_data = []
    index = []
    for tissue, ci in tissue_cis.items():
        tissue_data = ci[parameter]
        table_data.append(tissue_data)
        index.append(tissue)

    df = pd.DataFrame(table_data, index=index)

    return df
