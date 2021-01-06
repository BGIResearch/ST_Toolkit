# ConvertBinData
Convert bin expression matrix into bin1 data

## Platform & Environment
* python3.6.5+

## Prerequisties

| Package  | Version  | Description                                                |
| -------- | -------- | ---------------------------------------------------------- |
| pandas   | <=1.0.1  | handle dataframe                                           |
| numpy    | 1.19.1   |                                                            |
| scipy    | 1.5.2    | generate sparse matrix                                     |
| opencv   | 4.4.0    | generate h5ad format file to store gene expression matrix  |

## Run

### Usage
```
from ConverBinData import ConvertBinData

binSize = 10
convert = ConvertBinData()
mergedf = convert.ConvertData(partfile, genefile, binSize)
```

Required parameters:
* partfile: Downloaded bin gene expression file.
* genefile: The complete bin1 gene expression file.
* binSize: The binsize of your downloaded bin gene expression file. (Cannot larger than 50.)


