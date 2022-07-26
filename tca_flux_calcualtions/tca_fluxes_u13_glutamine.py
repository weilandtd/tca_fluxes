# -*- coding: utf-8 -*-
"""
This script contains all necessary instructions to fit the tissue TCA fluxes
from the provided U13-glutamine labeling data
"""


import pyximport
pyximport.install(setup_args={"script_args" : ["--verbose"]})

import sys
sys.path.append('../')

import multiprocessing
import json
import tqdm
import pandas as pd

from models.small_model import tca_model, IXM,  PARAMETERS

from tca_inference import InstatFluxFitter

"""
Computational parameters:

N_init  Number of initializations
N_CPU   Numper of CPUs used for the optimization (this should be choosen carefully)
"""

N_init = 10
N_CPU = 8

# Main Script
if __name__ == '__main__':

    # Load data sets (labeling and pools sizes)
    labeling_data = pd.read_csv('./../data/gln_labeling.csv')
    pool_size_data = pd.read_csv('./../data/pool_sizes.csv', )

    # Load initial guesses for the timescale
    tissue_time_scale = json.load(open('./../data/time_scale_estimates.json'))
    tissues = labeling_data.tissues.unique()


    # Initialize fitter class
    fitter = InstatFluxFitter(tca_model=tca_model,
                              metabolite_ix_hash_map=IXM,
                              parameters=PARAMETERS,)

    # Optimize each tissue
    for tissue in tissues:
        print("Running flux inference for {}".format(tissue))

        # Generate iterable input data for mulitprocessing
        fitting_inputs = fitter.make_input_data(tissue, pool_size_data, labeling_data, tissue_time_scale, N_init)

        # Run flux fitting using multiprocessing
        with multiprocessing.get_context('spawn').Pool(N_CPU) as pool:
            #TQDM for fancy progress bar
            results = list(tqdm.tqdm(pool.imap(fitter.fit_data, fitting_inputs), total=N_init))

        # Export results to a Table (pd.DataFrame)
        tissue_result = fitter.export_result(results, tissue=tissue)

        print("Postprocessing fitting results for {}".format(tissue))

        confidence_intervals = fitter.find_confidence_intervals(tissue_result, columns=['TCA','r'])