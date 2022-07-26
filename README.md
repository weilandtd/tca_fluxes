# Tools to determine TCA flux from in-stationary 13C infusions

This repository contains the tools and scripts to reproduce the TCA 
flux calculations presented in:

> Slow TCA flux implies low ATP production in tumors 
> Caroline R. Bartman, Yihui Shen, Won Dong Lee, Tara TeSlaa, 
> Connor S.R. Jankowski, Lin Wang, Lifeng Yang, Asael Roichman,
> Vrushank Bhatt, Taijin Lan, Zhixian Hu, Xi Xing, Wenyun Lu, 
> Jessie Yanxiang Guo, Joshua D. Rabinowitz
> doi: https://doi.org/10.1101/2021.10.04.463108

## Requirements
The following python packages are required to run the code: 

 - numpy
 - pandas
 - scipy
 - scikits.odes
 - matplotlib
 - cpython
 - tqdm
 - ipython (optional)

We recommend using [Anaconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
or [Miniconda](https://docs.conda.io/projects/conda/en/latest/glossary.html#miniconda-glossary)
to run the code, all required packages can be install by 
simply running: 
```bash
conda install --file requirements.txt
```

Alternatively you can use pip to install the required packages:

```bash
pip install -r requirements.txt
```
Note that when using pip installation you will need to install 
the [SUNDIALS](https://computing.llnl.gov/projects/sundials) 
solver suite as described 
[here](https://scikits-odes.readthedocs.io/en/stable/).

## Usage 
To use the scripts and tools provided it is first necessary to compile the models using `cython`:
```bash
cd path/to/repo/tca_fluxes/models
python setup.py build_ext --inplace
```
More details on how to use the full collision model and how to compile custom models can be found 
[here](https://github.com/weilandtd/tca_fluxes/tree/main/models). 

We then provide two scripts to perform the TCA flux inference using data from U13-Lactate 
(`tca_flux_calculations/tca_fluxees_u13_lactate.py`) 
and U13-Glutamine tracer (`tca_flux_calculations/tca_fluxees_u13_glutamine.py`)
Both these scripts follow a similar structure, first we import the required data which consists of 
three elements i) labeling data, ii) pool size data (tissue concentrations), and iii)
an estimate of the tissue timescale. 

```python
# Load data sets (labeling and pools sizes)
labeling_data = pd.read_csv('./../data/lac_labeling.csv')
pool_size_data = pd.read_csv('./../data/pool_sizes.csv', )

# Load initial guesses for the timescale
tissue_time_scale = json.load(open('./../data/time_scale_estimates.json'))
```

Second we initialize the parameter fitting class which requires 
i) a model function e.g. ii) a metabolite hash map, and 
iii) a list of the model parameter names these items can be imported from the respective model 
module: `from small_model import tca_model, IXM, PARAMETERS`

```python
fitter = InstatFluxFitter(tca_model=tca_model,
                          metabolite_ix_hash_map=IXM,
                          parameters=PARAMETERS,)
```

Third we run the actual parameter fitting for each tissue:

```python
# Generate iterable input data for mulitprocessing
fitting_inputs = fitter.make_input_data(tissue, pool_size_data, labeling_data, 
                                        tissue_time_scale, N_init)

# Run flux fitting using multiprocessing
with multiprocessing.get_context('spawn').Pool(N_CPU) as pool:
    #TQDM for fancy progress bar
    results = list(tqdm.tqdm(pool.imap(fitter.fit_data, fitting_inputs), total=N_init))

# Export results to a Table (pd.DataFrame)
tissue_result = fitter.export_result(results, tissue=tissue)
```

Finally, we determine the confidence intervals for the desired parameters here for the TCA turning flux`TCA`,
the malic enzyme flux`ME` and the PDH contribution to the TCA `r`:
```python
tissue_ci = fitter.find_confidence_intervals(tissue_result, columns=['TCA','ME','r'])

confidence_intervals[tissue] = tissue_ci
```

Using the `ci_table`function confidence intervals can be cast per parameter to a table:
```python
tca_ci = ci_table(confidence_intervals, 'TCA')
```