#! /usr/bin/env python3
# _*_ coding: utf-8 _*_

#@Author:  LiuXing liuxing2@genomics.cn
#@Date:  2021-01-06 15:45:38
#@Last Modified by: LiuXing
#@Last Modified time: 2021-01-06 15:45:38

import os, sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import logging
from optparse import OptionParser
import time

LOG_FORMAT="%(asctime)s %(levelname)s - %(message)s"
logging.basicConfig(level=logging.DEBUG, format=LOG_FORMAT)
TIME_FORMAT = "%y-%m-%d %H:%M:%S"

def main():
    actions = ["loom2h5ad", "markerGene2h5ad"]
    """
    %prog action [options]
    """
    parser = OptionParser(main.__doc__)
    parser.add_option("-i", "--in", action = "store", type = "str", dest = "inFile", help = "input file path.")
    parser.add_option("--h5ad", action = "store", type = "str", dest = "h5ad", help = "h5ad format output file path.")

    opts, args = parser.parse_args()

    if (len(args) < 1):
        print("please give the correct action, choose one from: " + ",".join(actions))
        sys.exit(not parser.print_help())
    if(opts.inFile == None):
        sys.exit(not parser.print_help())

    action = args[0].upper()
    if (action == "LOOM2H5AD"):
        loom2h5ad(opts.inFile, opts.h5ad)
    elif (action == "MARKERGENE2H5AD"):
        markerGene2h5ad(opts.inFile, opts.h5ad)
    else:
         print("please give the correct action, choose one from: " + ",".join(actions)) 
   
def loom2h5ad(loomFile, h5adFile):
    adata = sc.read_loom(loomFile, sparse=True, X_name='spliced', obs_names = 'cell', var_names='geneID', dtype = 'float32')
    adata.var.set_index(['Gene'], inplace=True)
    adata.obs['CellID'] = adata.obs['CellID'].apply(lambda x: x[1:].replace('_', '-'))
    adata.obs.set_index(['CellID'], inplace=True)
    adata.obs.rename(columns = {'seurat_clusters' : 'louvain', 'nFeature_RNA': 'n_genes', 'nCount_RNA': 'n_counts', 'coor_x': 'x', 'coor_y': 'y'}, inplace=True)
    adata.obsm['X_pca'] = adata.obsm.pop("pca_cell_embeddings")
    adata.obsm['X_umap'] = adata.obsm.pop("umap_cell_embeddings")
    adata.varm['PCs'] = adata.varm.pop("pca_feature_loadings")
    del adata.layers['scale_data']
    rawdata = anndata.AnnData(X = adata.X, obs = adata.obs[['coor_x','coor_y']], var = adata.var[['Selected']])
    adata.raw = rawdata
    adata.layers['raw_data'] = adata.X
    adata.X = adata.layers['norm_data']
    del adata.layers['norm_data']
    adata.write(h5adFile)

def markerGene2h5ad(csvFile, h5adFile):
    df = pd.read_csv(csvFile, index_col = 0, sep = "\s+")
    clusters = df['cluster'].nunique()
    dataList = [list() for i in range(clusters)]

    markerGenes = {'names': pd.DataFrame(), 'pvals': pd.DataFrame(), 'pvals_adj': pd.DataFrame(), 'pts': pd.DataFrame(), 
    'pts_rest': pd.DataFrame(), 'logfoldchanges': pd.DataFrame(), 'pts_dif': pd.DataFrame()}
    
    for i in range(clusters):
        dfC = df[df['cluster'] == i].reset_index()
        
        markerGenes['names'] = pd.concat([markerGenes['names'], dfC['gene'].to_frame()], axis=1)
        markerGenes['pvals'] = pd.concat([markerGenes['pvals'], dfC['p_val'].to_frame()], axis=1)
        markerGenes['pvals_adj'] = pd.concat([markerGenes['pvals_adj'], dfC['p_val_adj'].to_frame()], axis=1)
        markerGenes['pts'] = pd.concat([markerGenes['pts'], dfC['pct.1'].to_frame()], axis=1)
        markerGenes['pts_rest'] = pd.concat([markerGenes['pts_rest'], dfC['pct.2'].to_frame()], axis=1)
        markerGenes['logfoldchanges'] = pd.concat([markerGenes['logfoldchanges'], dfC['avg_logFC'].to_frame()], axis=1)
        markerGenes['pts_dif'] = pd.concat([markerGenes['pts_dif'], pd.DataFrame({'pts_dif': dfC['pct.1'] > dfC['pct.2']}).astype(int)], axis=1)
    
    dtypes = []*clusters
    for i in range(clusters):
        dtypes.append((str(i), 'O'))

    keys = ['names', 'pvals', 'pvals_adj', 'pts', 'pts_rest', 'logfoldchanges', 'pts_dif']
    adata = anndata.read(h5adFile)
    adata.uns['rank_genes_groups'] = dict()
    #print (markerGenes['names'])
    for key in keys:
        values = []
        for index, row in markerGenes[key].iterrows():
                values.append(tuple(row))
        adata.uns['rank_genes_groups'][key] = np.array(values, dtype = dtypes)
    adata.write(h5adFile) 

if __name__ == "__main__":
    main()