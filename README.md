# ST_Toolkit
tools to handle spatial transcriptome data

## compile
## Platform & Environment
* python3.6.5+

## Prerequisties

| Package  | Version  | Description                                                |
| -------- | -------- | ---------------------------------------------------------- |
| pandas   | 1.1.1.   | handle dataframe                                           |
| numpy    | 1.19.1   |                                                            |
| scipy    | 1.5.2    | generate sparse matrix                                     |
| anndata  | 0.7.4    | generate h5ad format file to store gene expression matrix  |
| optparse | 0.1.1    | manage command options.                                    |
| scanpy   | 1.6.0.   | read h5ad format file and do cell cluster                  |

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
                        The bin size that to combine the dnbs. default=50
  -t THREAD, --thread=THREAD
                        number of thread that will be used to run this
                        program. default=2
```
Required parameters:
* -i filename. Input gene expression matrix file in tsv format
* -0 file or directory name. Output path.
