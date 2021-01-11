#!/usr/bin/env python3

# @Author: LiuXing liuxing2@genomics.cn 
# @Date: 2020-12-03 11:00:53 
# @Last Modified by: LiuXing
# @Last Modified time: 2020-12-03 14:27:11
 

import os, sys
import pandas as pd
import numpy as np
import scipy
from scipy import ndimage
import matplotlib as mpl
import matplotlib.pyplot as plt
import json
from scipy import sparse
import logging
import numexpr
from multiprocessing import Pool
from optparse import OptionParser
import time

LOG_FORMAT="%(asctime)s %(levelname)s - %(message)s"
logging.basicConfig(level=logging.DEBUG, format=LOG_FORMAT)
TIME_FORMAT = "%y-%m-%d %H:%M:%S"

def main():
    actions = ["tsv2h5ad", "cellCluster", "visualization", "convertBinData"]
    """
    %prog action [options]
    """
    parser = OptionParser(main.__doc__)
    parser.add_option("-i", "--in", action = "store", type = "str", dest = "inFile", help = "input gene expression matrix file path.")
    parser.add_option("-o", "--out", action = "store", type = "str", dest = "out", help = "output file or directory path.")
    parser.add_option("-m", "--mask", action = "store", type = "str", dest = "mask", help = "lasso bin gene expression matrix get from stereomic visualization system.")
    parser.add_option("-s", "--binSize", action = "store", type = "int", dest = "binSize", default = 50, help = "The bin size or max bin szie that to combine the dnbs. default=50")
    parser.add_option("-t", "--thread", action = "store", type = "int", dest = "thread", default = 2, help = "number of thread that will be used to run this program. default=2")
    parser.add_option("-w", "--progress", action = "store", type = "int", dest = "progress", default = 4, help = "number of progress that will be used to run this program, only useful for visulization. default=4")

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
        slideBin.bin_stat()
    elif (action == "CELLCLUSTER"):
        cellCluster = CellCluster(opts.inFile, opts.out)
        cellCluster.scanpyCluster()
    elif (action == "VISUALIZATION"):
        visualization = Visualization(opts.inFile, opts.out, opts.binSize, opts.progress)
        visualization.process()
    elif (action == "CONVERTBINDATA"):
        convertBinData = ConvertBinData()
        convertBinData.ConvertData(opts.mask, opts.inFile, opts.out, opts.binSize)
    else:
        raise Exception("invalide action", 3)

class SlideBin():
    def __init__(self, geneExpFile, outdir, binSize):
        self.geneExpFile = geneExpFile
        self.outdir = outdir
        self.binSize = binSize

        os.makedirs(self.outdir, exist_ok=True)

    def bin_stat(self):
        import anndata
        import scanpy as sc
        df = pd.read_csv(self.geneExpFile, sep="\t")
        if "MIDCounts" in df.columns:
            df.rename(columns={"MIDCounts": "UMICount"}, inplace=True)
        
        df['x'] = (df['x']/self.binSize).astype(np.uint32)*self.binSize
        df['y'] = (df['y']/self.binSize).astype(np.uint32)*self.binSize
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
        logging.info("anndata generate finished...")

        resultFile = os.path.join(self.outdir, "{0}x{0}_adata.h5ad".format(self.binSize))
        adata.write(resultFile)

class CellCluster():
    def __init__(self, geneExpFile, outFile):
        self.geneExpFile = geneExpFile
        self.outFile = outFile

    def scanpyCluster(self):
        import scanpy as sc

        if (self.geneExpFile.endswith(".h5ad")):
            adata = sc.read_h5ad(self.geneExpFile)
        else:
            raise Exception("invalide file format", 1)

        sc.pp.filter_cells(adata, min_genes=0)
        sc.pp.filter_genes(adata, min_cells=0)
        adata.var['mt'] = adata.var_names.str.startswith(('mt-', 'MT-'))
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

        adata.raw = adata
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        if adata.var.shape[0] < 2000:
            raise Exception("there are no enough genes for cluster.", 2)
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
        adata = adata[:, adata.var.highly_variable]
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, svd_solver="arpack")
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
        sc.tl.tsne(adata)
        sc.tl.umap(adata)
        sc.tl.leiden(adata)
        sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
        adata.write(self.outFile)

class Visualization():
    def __init__(self, geneExpFile, outdir, maxBinSize=1000, progress=4):
        self.geneExpFile = geneExpFile
        self.outdir = outdir
        self.maxBinSize = maxBinSize
        self.progress = progress
        self.geneDf = pd.read_csv(self.geneExpFile, sep="\t")
        if "MIDCounts" in self.geneDf.columns:
            self.geneDf.rename(columns={"MIDCounts": "values"}, inplace=True)
        else:
            self.geneDf.rename(columns = {"UMICount": "values"}, inplace=True)
        dtypes = {"x": np.uint32, "y": np.uint32, "geneID": object, "values": np.uint32}
        self.geneDf.astype(dtypes, copy=False)
        os.makedirs(self.outdir, exist_ok=True)
        self.geneOutdir = os.path.join(outdir, "gene_merge")
        self.dnbOutdir = os.path.join(outdir, "dnb_merge")
        os.makedirs(self.geneOutdir, exist_ok=True)
        os.makedirs(self.dnbOutdir, exist_ok=True)

    def _get_geneTable(self):
        geneTable = self.geneDf['values'].groupby(self.geneDf['geneID']).sum().to_frame()
        geneTable.reset_index(inplace = True)
        geneTable.to_pickle(os.path.join(self.geneOutdir, "gene_table.pickle"))

    def _get_dnbRange(self):
        dnb_range_file = os.path.join(self.dnbOutdir, "dnb_range.json")
        dnb_range_map = {"min_x": int(self.geneDf['x'].min()), "max_x": int(self.geneDf['x'].max()), "min_y": int(self.geneDf['y'].min()), "max_y": int(self.geneDf['y'].max())}
        with open(dnb_range_file, 'w') as jsonHandle:
            json.dump(dnb_range_map, jsonHandle)

    def _get_dnb_pickle(self, geneBinDf, binSize):
        logging.info("start dnb pickle generation of bin{0}".format(binSize))
        start = time.time()
        dnbDf = geneBinDf['values'].groupby([geneBinDf['x'], geneBinDf['y']]).sum().to_frame()
        dnbDf['gene_counts'] = geneBinDf['geneID'].groupby([geneBinDf['x'], geneBinDf['y']]).nunique()
        dnbDf.reset_index(inplace = True)
        binFile = os.path.join(self.dnbOutdir, "merge_dnb_bin{0}.pickle".format(binSize))
        dnbDf.to_pickle(binFile)
        if (binSize==200):
            self._get_fig(dnbDf, 200)
        logging.info("finish dnb pickle generation of bin{0}".format(binSize))
        end = time.time()
        logging.info("time used for bin{0} generation: {1}".format(binSize, end-start))

    def _get_fig(self, dnbDf, binSize):
        cmap = mpl.colors.ListedColormap(['#0C3383', '#0A88BA', '#F2D338', '#F28F38', '#D91E1E'])
        x_range = dnbDf['x'].max() - dnbDf['x'].min()
        y_range = dnbDf['y'].max() - dnbDf['y'].min()

        plt.figure(figsize=(1,y_range/x_range))
        plt.scatter(dnbDf['x'], dnbDf['y'], c=dnbDf['values'], s=1, cmap=cmap)
        plt.axis('off')
        plt.savefig(os.path.join(self.dnbOutdir, "bin{0}.png".format(binSize)))

    def _get_gene_pickle(self, geneBinDf, binSize):
        logging.info("start gene pickle generation of bin{0}".format(binSize))
        start = time.time()
        if (binSize == 1):
            groupDf = geneBinDf
        else:
            geneBinDf['x'] = geneBinDf['x'].map(lambda x: int(x/binSize)*binSize).astype(np.uint32)
            geneBinDf['y'] = geneBinDf['y'].map(lambda y: int(y/binSize)*binSize).astype(np.uint32)
            geneBinDf['geneID'] = geneBinDf['geneID'].astype(str)
            groupDf = geneBinDf['values'].groupby([geneBinDf['x'], geneBinDf['y'], geneBinDf['geneID']]).sum().reset_index()
            groupDf['geneID'] = groupDf['geneID'].astype('category')
        geneBinFile = os.path.join(self.geneOutdir, "merge_gene_bin{0}.pickle".format(binSize))
        groupDf.to_pickle(geneBinFile)
        logging.info("finish gene pickle generation of bin {0}".format(binSize))
        end = time.time()
        logging.info("time used for bin{0} generation: {1}".format(binSize, end-start))
        return groupDf
    
    def _get_pickle(self, geneBinDf, binSize):
        groupDf = self._get_gene_pickle(geneBinDf, binSize)
        self._get_dnb_pickle(groupDf, binSize)

    def process(self):
        binSizeList = [1, 2, 5, 8, 10, 15, 20, 30, 50, 80, 100, 120, 150, 180, 200, 250, 300, 400, 500, 600, 800, 1000]
        binSizeList = filter(lambda x: x<=self.maxBinSize, binSizeList)
        self._get_geneTable()
        self._get_dnbRange()
       
        pool = Pool(self.progress)
        for binSize in binSizeList:
            pool.apply_async(self._get_pickle, (self.geneDf, binSize, ))
        pool.close()
        pool.join()

class ConvertBinData():
    """
    # @ Author: Xiaoxuan Tang
    # @ Create Time: 2021-01-06 11:26:59
    # @ Modified by: xiaoxuan Tang
    # @ Modified time: 2021-01-06 14:26:05
    input: The lasso bin gene expression matrix; The complete gene expression matrix
    return: Binsize=1 gene expression matrix.
    """

    def __init__(self):
        self.typeColumn = {"geneID": 'str', "x": np.uint32, \
                            "y": np.uint32, "values": np.uint32, 'MIDCount':np.uint32, \
                            "MIDCounts":np.uint32, "UMICount": np.uint32}
        
    def __Dumpresult(self, mask, genedf):
        dst = np.where(mask > 0)
        dstx = dst[1]
        dsty = dst[0]
        tissue = pd.DataFrame()
        tissue['x'] = [ii + self.Xmin for ii in dstx]
        tissue['y'] = [ij + self.Ymin for ij in dsty]

        mergedf = pd.merge(genedf, tissue, how='inner', on=['x', 'y'])

        return mergedf

    def __CreateImg(self, df):
        bindf = pd.DataFrame()
        bindf['x'] = df['x'] - self.Xmin
        bindf['y'] = df['y'] - self.Ymin
        bindf['values'] = [255] * len(df)
        
        sparseMt = sparse.csr_matrix((bindf['values'].astype(np.uint8), (bindf['y'], bindf['x'])))
        img = sparseMt.toarray()
        return img

    def ConvertData(self, partfile, genefile, outFile, binSize):
        import cv2
        ### Initial data
        if binSize > 50:
            raise ValueError("Binsize could not larger than 50.")

        genedf = pd.read_csv(genefile, sep='\t', dtype=self.typeColumn)
        partdf = pd.read_csv(partfile, sep='\t', dtype=self.typeColumn)
        
        if len(genedf) < len(partdf):
            raise Warning("Inputs are not correct.")
        self.Xmin = genedf['x'].min()
        self.Ymin = genedf['y'].min()
        print("Processing data..")
        ### Create Mask
        part_img = self.__CreateImg(partdf)
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2*binSize, 2*binSize))
        mask = cv2.morphologyEx(part_img, cv2.MORPH_CLOSE, kernel)

        ### Dump results
        mergedf = self.__Dumpresult(mask, genedf)
        mergedf.to_csv(outFile, sep="\t")
        #return mergedf

if __name__ == "__main__":
    main()
