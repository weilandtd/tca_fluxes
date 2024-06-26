# Tools to determine TCA flux from non-stationary 13C infusions

This repository contains the tools and scripts to reproduce the TCA 
flux calculations presented in:

> Caroline R. Bartman, Daniel R. Weilandt, Yihui Shen, Won Dong Lee, Yujiao Han,
> Tara TeSlaa, Connor S.R. Jankowski, Laith Samarah, Noel Park, Maria Victoria da Silva,
> Maya Aleksandrova, Yetis Gultekin, Lin Wang, Lifeng Yang, Asael Roichman, Vrushank Bhatt,
> Taijin Lan,Zhixian Hu, Xi Xing, Wenyun Lu, Shawn Davidson, Matthew Vander Heiden, Daniel Herranz,
> Jessie Yanxiang Guo, Yibin Kang, Joshua D. Rabinowitz, **Slow TCA flux and ATP production in primary solid tumors**. *Nature* (2023).
> https://doi.org/10.1038/s41586-022-05661-6


We here implement non-stationary metabolic flux analysis using a simplified TCA model (see figure below). 
The code allows the user to infer the fluxes of reactions in the network depicted below from labeling data of
non-stationary isotope infusions. A full description of the methodology can be found 
[here](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05661-6/MediaObjects/41586_2022_5661_MOESM1_ESM.pdf).

![Model schematic](https://github.com/weilandtd/tca_fluxes/blob/main/_img.png)


## Download
To download the data and code from this repository use git clone and make sure you have [GIT-LFS](https://git-lfs.github.com/) installed:

```bash
git clone https://github.com/weilandtd/tca_fluxes
git lfs install 
git pull
```


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
conda config --add channels conda-forge
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
(`tca_flux_calculations/tca_fluxes_u13_lactate.py`) 
and U13-Glutamine tracer (`tca_flux_calculations/tca_fluxes_u13_glutamine.py`)
To run these scripts change the directory to `tca_flux_calculations` and open an `ipython` console. 
You can then run the scripts:
```bash
run example_tca_fluxees_u13_lactate.py
```
when the script starts successful you should see an output similar to below. 
Note: Some of the initial guesses may result in $[CVODE ERROR] ... the corrector convergence test failed repeatedly ...$. Since we use random sampling over a wide range of parameters this cannot be avoided but does not mean that the algorism is not working.
```
Running flux inference for diaphragm
100%|██████████████████████████████| 1000/1000 [13:42<00:00,  1.22it/s]
Postprocessing fitting results for diaphragm
...
```
Note that the scripts `tca_fluxes_u13_lactate.py` and `tca_fluxes_u13_glutamine.py` are configured 
to reproduce the paper results and are quite computationally expensive (See code for details).  

## Additional information

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
