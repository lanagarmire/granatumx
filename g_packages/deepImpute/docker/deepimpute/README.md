# DeepImpute: an accurate and efficient deep learning method for single-cell RNA-seq data imputation 
Arisdakessian, C et al. "DeepImpute: an accurate and efficient deep learning method for single-cell RNA-seq data imputation", BioRxiv (2018): XXXX\
DeepImpute has been implemented in python2 and python3. The recommanded version is python3.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Installing

To install DeepImpute, you only need to download the git repository at https://github.com/lanagarmire/DeepImpute and install it using pip:

```
>> git clone https://github.com/lanagarmire/DeepImpute
>> cd DeepImpute
>> pip install -e . 
```

### Usage 

DeepImpute can be used either in commandline or as a python package.

Command line: 
```

usage: deepImpute.py [-h] [-o outputFile] [--ncores NCORES]
                     [--cell-axis {rows,columns}] [--limit LIM]
                     [--subset SUBSET]
                     inputFile

scRNA-seq data imputation using DeepImpute.

positional arguments:
  inputFile             Path to input data.

optional arguments:
  -h, --help            show this help message and exit
  -o outputFile         Path to output data counts. Default: ./
  --cores NCORES       Number of cores. Default: 5
  --cell-axis {rows,columns}
                        Cell dimension in the matrix. Default: rows
  --limit LIM             Genes to impute (e.g. first 2000 genes). Default: auto
  --subset SUBSET       Cell subset to speed up training. Either a ratio
                        (0<x<1) or a cell number (int). Default: 1 (all)
```

Python package:
```
from deepimpute.deepImpute import deepImpute

data = pd.read_csv('../examples/test.csv',index_col=0) # dimension = (cells x genes)
imputed = deepImpute(data,NN_lim='auto',n_cores=10,cell_subset=1)
```

A more detailed usage of deepImpute's functionnalities is available in the iPython Notebook notebook\_example.ipynb


### Running the tests

Each file has been validated using a unittest script. They are all available in the test folder.\n\n 
To run all the tests at once, you can also use the makefile by running ```make test```.

## License

