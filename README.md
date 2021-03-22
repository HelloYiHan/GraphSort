# GraphSort
GraphSort is a novel tool for in silico estimation of immune and tissue cells compositions. Gene interactive relation information in the Kyoto Encyclopedia of Genes and Genomes (KEGG) is extracted by GraphSort with graph convolutional network (GCN) and it helped improve the accuracy of estimation. The signature matrix is no longer needed anymore for GraphSort thanks to the GCN. Compared with four mainstream and up-to-date deconvolution software, GraphSort achieves a higher Pearson correlation and lower absolute error of estimated cell fractions with respect to the ground truth.

## Brief Introduction
The workflow of GraphSort:
![alt text](https://github.com/HelloYiHan/GraphSort/blob/master/FigInGitHub.png?raw=true)

## Run Online
GraphSort could be installed and run on [Google Colab online](https://colab.research.google.com/drive/1n8IYkP8SeSOhyrYDdcJYJ6EtMwK7NWBn?usp=sharing)

Please read the manual below, especially the Formatting Requirements of Input Mixture File part before running the program.

**Note**: Sometimes the required python package torch-sparse==0.4.0 can not be successfully installed since its developer no longer maintains this old version any more. We provide the docker image of our GraphSort instead.

## Docker
First, download the  docker image.

Second, load the docker image.
```shell
docker load --input ./GraphSortV101.tar
```

Third, test an example dataset. **Remeber** to replace the YourWorkingDirectory(1 places) with your path.

```shell
docker run -it -v YourWorkingDirectory:/local graphsort:v1.0.1 python /workspace/GraphSort-master/run_graphsort.py --input /workspace/GraphSort-master/example_data_gse107011.txt --type rnaseq --output /local/test_107011_out.txt
```

Fourth, run your datasets. **Remeber** to replace the YourData(2 places) with your file name and the YourWorkingDirectory(1 places) with your path.

```shell
docker run -it -v YourWorkingDirectory:/local graphsort:v1.0.1 python /workspace/GraphSort-master/run_graphsort.py --input /local/YourData.txt --type rnaseq/microarray/pancreatic --output /local/YourData_out.txt
```

## Local Installation

### PyTorch and Extension Libraries
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
Preprocessing is implemented in R and the following packages are needed: KEGGgraph, KEGG.db, testit, SparseM, graph, funr, edgeR, sva, and preprocessCore.

### Dependent Datasets
Training datasets are needed to remove the batch effect of input and training datasets.

1. [Download](https://drive.google.com/file/d/18DoMwpMa8PFajx_Q-gXWMrWOsEYBkXtE/view?usp=sharing) the RNA-Seq training dataset with 1k samples (simulated).
2. [Download](https://drive.google.com/file/d/19K_qwpuI5eHPr1l5He3p0vY8wt3Uq7sH/view?usp=sharing) the RNA-Seq training dataset with 2k samples (simulated). It uses more memory than the 1k training dataset.
3. [Download](https://drive.google.com/file/d/1NnlqQbPd2xC7lHSZaeVWTz6MuBhI_RLm/view?usp=sharing) the Microarray training dataset.
4. [Download](https://drive.google.com/file/d/1erBWx8JNFzfePlIXriulrghxaaKAHxl7/view?usp=sharing) the pancreatic islet training dataset.
5. Uncompress these datasets and put them in the GraphSort directory, e.g.
```
unzip rem_bat_eff_dat_n1000.zip
mv rem_bat_eff_dat_n1000.txt ./GraphSort
```

## Usage
```
python run_graphsort.py [arguments]

Required arguments:
--input        -i      The expression file to be analyzed.
--type         -t      The type of the input file: 'rnaseq' OR 'microarray' for immune cells; 'pancreatic' for endocrine cells in the human pancreatic islets.

Optional arguments:
--output       -o      Filename of the estimated output. Default: graphsort_out.txt
--batch_size   -b      Batch size for computation. Default: 10
--device       -d      Computation device: gpu OR cpu. Default: cpu
--size         -s      Data size in preprocessing(only matters for RNA-Seq data): 1k OR 2k. Default: 1k
```
## Formatting Requirements of Input Mixture File 
The `format` of the input `RNA-Seq` data should be:
* `Tab-delimited` with no double quotations and no missing entries.
* The first row consists of sample names.
* The first column consists of `Ensembl Gene IDs`.
* The values are not log-transformed.
* The values should be raw counts.

The `format` of the input `Microarray` data should be:
* `Tab-delimited` with no double quotations and no missing entries.
* The first row consists of sample names.
* The first column consists of `Gene Symbols`.
* The values are not log-transformed.

The `format` of  `human pancreatic islets` data should be:
* `Tab-delimited` with no double quotations and no missing entries.
* The first row should be sample names.
* The first column should be the `Gene Symbols`.
* The values are not log-transformed and raw RNA-Seq counts.

## Output
To estimate immune cells with RNA-Seq data, GraphSort could estimate fractions of 7 immune cell types: B cells, CD4 T cells, CD8 T cells, Monocytes, Basophils, Dendritic cells, and NK cells.

To estimate immune cells with Microarray data, GraphSort could estimate fractions of 8 immune cell types: Memory B cells, Naive B cells, Plasma cells, CD4 T cells, CD8 T cells, Monocytes, Dendritic cells, and NK cells.

To estimate endocrine cells in human pancreatic islets, GraphSort could estimate fractions of 4 endocrine cell types: alpha cells, beta cells, gamma cells, and delta cells. The exocrine ductal and acinar cells are also estimated.

## Example
A detailed example about downloading datasets from GEO, preprocessing the input file (R code), uploading the input file to Google Colab, and running GraphSort is in the GitHub directory (detailed example.txt).

Example datasets [GSE107011](https://drive.google.com/file/d/1nqR8C_x8uX8J6ifZqmfnF29cLTvRkhkQ/view?usp=sharing), [SDY67](https://drive.google.com/file/d/1tTb4dSL_xj-eEoaRuu27faOX07gRBNNX/view?usp=sharing), [GSE65133](https://drive.google.com/file/d/1j9taMZvFkD8mTs0hLjRLn0PijLyIhI0z/view?usp=sharing), [GSE106898](https://drive.google.com/file/d/1eLcqgCdDG7aVoa-cNpSORG5wkYPPgyuz/view?usp=sharing), [E-MTAB5060](https://drive.google.com/file/d/1tS6RdUj93iKTzYfjAq7Eq0Y-SBQXjxcB/view?usp=sharing), [GSE59654](https://drive.google.com/file/d/1T4CROU-iJX3mZfKAhNsxsO4sxuFhcwnJ/view?usp=sharing), and [GSE60424](https://drive.google.com/file/d/1_TEzuku9ePRiKFpHUzbsQVK8xEX5Er1A/view?usp=sharing) could be downloaded here.

Below is the running command of GraphSort:

```
python ./GraphSort/run_graphsort.py --input example_data_gse107011.txt --type rnaseq --output gse107011_out.txt
```

## Citation
Please cite our paper if you use GraphSort in your work:
> Han, Y. et al. (2020) GraphSort: a geometric deep learning algorithm for in silico dissecting cell compositions in bulk expression data. Manuscript submitted for publication.
