import pyximport

pyximport.install(setup_args={"script_args": ["--verbose"]})

import sys

sys.path.append('../')

import multiprocessing
import json
import tqdm
import pandas as pd
import numpy as np

from models.small_model import tca_model, IXM, PARAMETERS

from tca_inference import InstatFluxFitter, ci_table

# Main Script
if __name__ == '__main__':
    fitting_results = pd.read_csv('lactate_parameter_estimates_2000.csv')

    # Select a parameters set
    single_fitting_result = fitting_results.loc[0]

    # Initialize fitter class
    fitter = InstatFluxFitter(tca_model=tca_model,
                              metabolite_ix_hash_map=IXM,
                              parameters=PARAMETERS, )

    parameters, pool_sizes = fitter.get_pools_and_params(single_fitting_result)

    # Time in minutes and number of timepints
    t_max = 100
    N = 100
    t = np.linspace(0, t_max, N)

    sol = fitter.solve_odes(t, parameters, pool_sizes)
    fitter.export_ode_solution(sol, folder='output', name=single_fitting_result['tissue'])
