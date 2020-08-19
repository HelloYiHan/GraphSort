# GraphSort

GraphSort is a novel tool for in silico estimation of immune cell compositions in blood. Gene interactive relation information is extracted by GraphSort with graph convolutional network (GCN) and helped improve the accuracy of estimation. Signature matrix is no longer needed any more for GraphSort thanks to the GCN.

## Brief introduction

work flow

## Installation

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


## Usage
```
python run_graphsort.py

Arguments
--input, -i  RNA-Seq or Microarray expression data of samples to be estimated.
```

## Examples
