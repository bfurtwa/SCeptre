# SCeptre
**(Single Cell proteomics readout of expression)**

SCeptre is a python package that extends the functionalities of Scanpy to analyze multiplexed single-cell proteomics data.

## Installation from Github

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

## Usage
### Replication of Schoof et al.

Usage is exemplified in the notebooks for the analysis from the mansucript "Quantitative Single-Cell Proteomics as a Tool to Characterize Cellular Hierarchies".

The analysis can be replicated using the provided conda environment:

```
conda env create -f Schoof_et_al/code/environment.yml
conda activate sceptre
pip install Schoof_et_al/code/sceptre-0.1-py3-none-any.whl
```

The required data can be downloaded from
http://proteomecentral.proteomexchange.org
using the dataset identifier PXD020586

Find the notebooks in the subdirectory `Schoof_et_al/code`, place the required data in `Schoof_et_al/data`, and create the folder `Schoof_et_al/results/tmp/`.


The following notebooks process the different datasets:

Notebook | Description
------------ | -------------
300ms.ipynb | Sceptre analysis of the 'medium' dataset.
500ms.ipynb | SCeptre analysis of the 'high' dataset.
bulk.ipynb | SCeptre analysis of the 'bulk' dataset.
integrated.ipynb | SCeptre analysis of the 'integrated' dataset.

### Functions and example worklow
Each function has its own docstring explaining the function in depth. A typical workflow makes usage of the following steps:

#### Create meta data
To create the meta data for each cell, as done in Schoof et al., from a collection of tables describing the experimental design
and layouts of the 384-well plates, the following function is used. For details on the required tables, have a look at `Schoof_et_al/data/500ms`.
```
import sceptre as spt
spt.create_meta_data(input_dir="../data/500ms/", output_dir=res_dir)
```
Alternatively, the meta data table can be created by the user. It requires the columns `File ID` and `Channel` to map the meta data to each cell.

#### Load dataset
To load the dataset into python, to following function is used. To this end, only output tables from Proteome Discoverer are supported.
```
dataset = spt.load_dataset(proteins = "../data/500ms/500ms_Proteins.txt",
                           psms = "../data/500ms/500ms_PSMs.txt",
                           msms = "../data/500ms/500ms_MSMSSpectrumInfo.txt",
                           files = "../data/500ms/500ms_InputFiles.txt",
                           meta = res_dir + "meta.txt")
```

#### LC-MS QC
The dataset object can be used to quality control each LC-MS run with the follwing functions.
```
spt.plot_psms_msms(dataset)
spt.plot_avg_sn(dataset)
spt.plot_set_overview(dataset)

s_c_channels = ['127N', '128N', '128C', '129N', '129C', '130N', '130C',
                '131N','131C', '132N', '132C', '133N', '133C', '134N']
spt.print_ms_stats(dataset, s_c_channels=s_c_channels)

spt.plot_interference(dataset)
```

#### Load dataset into Scanpy
Subsequently the dataset object is used to create a scanpy adata object.
```
adata = spt.dataset_to_scanpy(dataset)
```

#### Filtering and normalization
Non-single cell channels have to be removed.
```
adata = adata[adata.obs['Channel'] != '126'].copy()
adata = adata[adata.obs['Channel'] != '127C'].copy()
```
Then the dataset can be normalized.
```
spt.normalize(adata)
```

#### Cell QC
The follwing functions are used to filter out outlier cells.
```
spt.calculate_cell_filter(adata)
spt.apply_cell_filter(adata)
```

#### Batch QC
To detect potential systematic bias introduced by the sample preparation or measurement the following functions can be used.
```
spt.plot_batch_qc(adata)
spt.plot_plate_qc(adata)
```

#### Scanpy analysis
The adata object can now be used to perform a scanpy analysis. 