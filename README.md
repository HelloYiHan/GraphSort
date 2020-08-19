# GraphSort
GraphSort is a novel tool for in silico estimation of immune cell compositions in blood. Gene interactive relation information is extracted by GraphSort with graph convolutional network (GCN) and helped improve the accuracy of estimation. Signature matrix is no longer needed any more for GraphSort thanks to the GCN.

## Brief introduction

work flow figure

## Installation
Setup a conda environment is recommended.
### PyTorch and extension libraries
GraphSort was developed with PyTorch 1.2.0 and [PyTorch Geometric](https://github.com/rusty1s/pytorch_geometric) 1.3.2.
```
pip install torch==1.2.0 torchvision==0.4.0
pip install torch-scatter==1.3.1
pip install torch-sparse==0.4.0
pip install torch-cluster==1.4.4
pip install torch-spline-conv==1.1.0
pip install torch-geometric==1.3.2
```
### Standard R Package
Preprocessing is implemented in R and the following packages are needed: KEGGgraph, KEGG.db, testit, SparseM, graph, funr, edgeR, sva, preprocessCore.

### Denendent datasets
Training datasets are needed to remove batch effect of input and training datasets.

1. [Download](https://drive.google.com/file/d/18DoMwpMa8PFajx_Q-gXWMrWOsEYBkXtE/view?usp=sharing) the RNA-Seq training dataset with 1k samples (simulated).
2. [Download](https://drive.google.com/file/d/19K_qwpuI5eHPr1l5He3p0vY8wt3Uq7sH/view?usp=sharing) the RNA-Seq training dataset with 2k samples (simulated). It use more memory than the 1k training dataset.
3. [Download](https://drive.google.com/file/d/1NnlqQbPd2xC7lHSZaeVWTz6MuBhI_RLm/view?usp=sharing) the Microarray training dataset.


4. Uncompress these datasets and  put them in the GraphSort directory.
e.g.
```
unzip rem_bat_eff_dat_n1000.zip
mv rem_bat_eff_dat_n1000.txt ./GraphSort
```
## Installation Online
Install in Google Colab online.


## Usage
```
python run_graphsort.py [arguments]

Required arguments:
--input        -i      Expression data (RNA-Seq or microarray) of samples to be estimated
--type         -t      Expression data type: rnaseq OR microarray

Optional auguments:
--output       -o      File name of the estimation output. Default: graphsort_out.txt
--batch_size   -b      Batch size for computation. Default: 10
--device       -d      Computation device: gpu OR cpu. Default: cpu
--size         -s      Data size in preprocessing(only matters for RNA-Seq data): 1k OR 2k. Default: 1k
```

## Examples
```
python ./GraphSort/run_graphsort.py --input example_data_gse107011.txt --type rnaseq --output gse107011_out.txt
```

## Cite
Please cite our paper if you use this code in your own work:
```
Yi Han, et al. (2020). GraphSort: a geometric deep learning algorithm for in silico dissecting cell compositions in bulk expression data. In preparation.
```
