# ST_Toolkit
tools to handle spatial transcriptome data

## compile
## Platform & Environment
* python3.6.5+

## Prerequisties

| Package  | Version  | Description                                                |
| -------- | -------- | ---------------------------------------------------------- |
| pandas   | <=1.0.1  | handle dataframe                                           |
| numpy    | 1.19.1   |                                                            |
| scipy    | 1.5.2    | generate sparse matrix                                     |
| anndata  | 0.7.4    | generate h5ad format file to store gene expression matrix  |
| optparse | 0.1.1    | manage command options.                                    |
| scanpy   | 1.6.0    | read h5ad format file and do cell cluster                  |
| numexpr  | 2.7.1    | accellerate operation on numpy array                       |

## Run

### Usage
```
python3 ST_Handle_Exp.py 

Usage: ST_Handle_Exp.py action [options]

Options:
  -h, --help            show this help message and exit
  -i INFILE, --in=INFILE
                        input gene expression matrix file path.
  -o OUT, --out=OUT     output file or directory path.
  -s BINSIZE, --binSize=BINSIZE
                        The bin size or max bin szie that to combine the dnbs.
                        default=50
  -t THREAD, --thread=THREAD
                        number of thread that will be used to run this
                        program. default=2
  -w PROGRESS, --progress=PROGRESS
                        number of progress that will be used to run this
                        program, only useful for visulization. default=4
```
Required parameters:
* -i INFILE, --in=INFILE input gene expression matrix file path.
* -o OUT, --out=OUT  output file or directory path.

Optional parameters:
* -s BINSIZE, --binSize=BINSIZE The bin size or max bin szie that to combine the dnbs.default=50.
* -t THREAD, --thread=THREAD number of thread that will be used to run this program. default=2.
* -w PROGRESS, --progress=PROGRESS number of progress that will be used to run this program, only useful for visulization. default=4.

#### action-tsv2h5ad
```
description:
  merge number of binSize*binSize dnb data to one spot and save new spatial gene expression matrix to h5ad format file. 
e.g:
  python3 ST_Handle_Exp.py tsv2h5ad -i merge_GetExp_gene.txt -o merge_GetExp_gene.h5ad -s 50
```   
#### action-visualization
```
discription:
  generate pickle format file with spatial gene expression matrix that has been grouped by different binSize.
e.g:
  python3 ST_Handle_Exp.py visualization -i merge_GetExp_gene.txt -o visualization -s 1000
```
