#!/usr/bin/env python3

# @Author: LiuXing liuxing2@genomics.cn 
# @Date: 2020-12-03 11:00:53 
# @Last Modified by: LiuXing
# @Last Modified time: 2020-12-03 14:27:11
 

import os, sys
import pandas as pd
import numpy as np
import anndata
import scanpy as sc
from scipy import sparse
import logging
import numexpr
from optparse import OptionParser

LOG_FORMAT="%(asctime)s %(levelname)s - %(message)s"
logging.basicConfig(level=logging.DEBUG, format=LOG_FORMAT)
TIME_FORMAT = "%y-%m-%d %H:%M:%S"

def main():
    actions = ["tsv2h5ad"]
    """
    %prog action [options]
    """
    parser = OptionParser(main.__doc__)
    parser.add_option("-i", "--in", action = "store", type = "str", dest = "inFile", help = "input gene expression matrix file path.")
    parser.add_option("-o", "--out", action = "store", type = "str", dest = "out", help = "output file or directory path.")
    parser.add_option("-s", "--binSize", action = "store", type = "int", dest = "binSize", default = 50, help = "The bin size that to combine the dnbs. default=50")
    parser.add_option("-t", "--thread", action = "store", type = "int", dest = "thread", default = 2, help = "number of thread that will be used to run this program. default=2")

    opts, args = parser.parse_args()

    if (len(args) < 1):
        print("please give the correct action, choose one from: " + ",".join(actions))
        sys.exit(not parser.print_help())
    if(opts.inFile == None):
        sys.exit(not parser.print_help())

    numexpr.set_num_threads(opts.thread)

    action = args[0].upper()
    if (action == "TSV2H5AD"):
        slideBin = SlideBin(opts.inFile, opts.out, opts.binSize)
        slideBin.BinStat() 

class SlideBin():
    def __init__(self, geneExpFile, outdir, binSize):
        self.geneExpFile = geneExpFile
        self.outdir = outdir
        self.binSize = binSize

        os.makedirs(self.outdir, exist_ok=True)

    def BinStat(self):
        df = pd.read_csv(self.geneExpFile, sep="\t")
        if "MIDCounts" in df.columns:
            df.rename(columns={"MIDCounts": "UMICount"}, inplace=True)
        
        df['x'] = (df['x']/self.binSize).astype(np.uint32)
        df['y'] = (df['y']/self.binSize).astype(np.uint32)
        df['cell'] = df['x'].astype(str) + "-" + df['y'].astype(str)
        bindf = df['UMICount'].groupby([df['cell'], df['geneID']]).sum()
        cells = set(x[0] for x in bindf.index)
        genes = set(x[1] for x in bindf.index)
        cellsdic = dict(zip(cells, range(0, len(cells))))
        genesdic = dict(zip(genes, range(0, len(genes))))
        rows = [cellsdic[x[0]] for x in bindf.index]
        cols = [genesdic[x[1]] for x in bindf.index]
        expMtx = sparse.csr_matrix((bindf.values, (rows, cols)))

        obs = pd.DataFrame(index = cells)
        var = pd.DataFrame(index = genes)
        adata = anndata.AnnData(X = expMtx, obs = obs, var = var)
        logging.info("ammdata generate finished...")

        resultFile = os.path.join(self.outdir, "{0}x{0}_adata.h5ad".format(self.binSize))
        adata.write(resultFile)

if __name__ == "__main__":
    main()
