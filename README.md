# SCeptre
**(Single Cell proteomics readout of expression)**

SCeptre is a python package that extends the functionalities of Scanpy to analyze multiplexed single-cell proteomics data.

### Installation from Github

Tested on Ubuntu 20.04.1 LTS.  
It's recommended to work in a dedicated conda environment. E.g:

```
conda create -n sceptre python=3.7
conda activate sceptre
```

Clone the repository and `cd` into its root directory. Then:

```
pip install .
```

### Usage

Usage is exemplified in the notebooks for the analysis from the mansucript "Quantitative Single-Cell Proteomics as a Tool to Characterize Cellular Hierarchies".

The analysis can be replicated using the provided conda environment:

```
conda env create -f Schoof_et_al/code/environment.yml`
conda activate sceptre
pip install Schoof_et_al/code/sceptre-0.1-py3-none-any.whl
```

Find the notebooks in the subdirectory `Schoof_et_al/code` and the required data in `Schoof_et_al/data`.

The following notebooks process the different datasets:

Notebook | Description
------------ | -------------
300ms.ipynb | Sceptre analysis of the 'medium' dataset.
500ms.ipynb | SCeptre analysis of the 'high' dataset.
bulk.ipynb | SCeptre analysis of the 'bulk' dataset.
integrated.ipynb | SCeptre analysis of the 'integrated' dataset.

