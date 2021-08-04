#!/usr/bin/env python3

# @Author: LiuXing liuxing2@genomics.cn 
# @Date: 2020-12-03 11:00:53 
# @Last Modified by: LiuXing
# @Last Modified time: 2020-12-03 14:27:11
 

import os, sys
import pandas as pd
import numpy as np
from pandas.core.reshape.merge import merge
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
import h5py
import math
import csv

LOG_FORMAT="%(asctime)s %(levelname)s - %(message)s"
logging.basicConfig(level=logging.DEBUG, format=LOG_FORMAT)
TIME_FORMAT = "%y-%m-%d %H:%M:%S"

def main():
    actions = ["tsv2h5ad", "dnbMatting", "cellCluster", "visualization", "convertBinData", "saturationPlot", "mergeGem"]
    """
    %prog action [options]
    """
    parser = OptionParser(main.__doc__)
    parser.add_option("-i", "--in", action = "store", type = "str", dest = "inFile", help = "input gene expression matrix file path.")
    parser.add_option("-o", "--out", action = "store", type = "str", dest = "out", help = "output file or directory path.")
    parser.add_option("-s", "--binSize", action = "store", type = "int", dest = "binSize", default = 50, help = "The bin size or max bin szie that to combine the dnbs. default=50")
    parser.add_option("-m", "--mask" , action = "store", type = "str", dest = "mask", help = "input gene expression matrix file generated from lesso tool of the stereomap system.")
    parser.add_option("--geneNumber", action = "store", type = "int", dest = "geneNumber", default = 2000, help = "number of genes will be used to cluster bins.")
    parser.add_option("-t", "--thread", action = "store", type = "int", dest = "thread", default = 8, help = "number of thread that will be used to run this program. default=2")
    parser.add_option("-w", "--progress", action = "store", type = "int", dest = "progress", default = 1, help = "number of progress that will be used to run this program, only useful for visulization. default=4")

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
        slideBin.process()
    elif (action == "CELLCLUSTER"):
        cellCluster = CellCluster(opts.inFile, opts.out, opts.binSize)
        cellCluster.scanpyCluster()
    elif (action == "DNBMATTING"):
        dnbMatting = DnbMatting(opts.inFile, opts.out, opts.binSize)
        dnbMatting.dnb_matting()
    elif (action == "VISUALIZATION"):
        visualization = Visualization(opts.inFile, opts.out, opts.binSize, opts.progress)
        visualization.process()
    elif (action == "CONVERTBINDATA"):
        convertBinData = ConvertBinData(opts.mask, opts.inFile, opts.out, opts.binSize)
        convertBinData.ConvertData()        
    elif (action == "SATURATIONPLOT"):
        saturationPlot = SaturationPlot(opts.inFile, opts.out)
        saturationPlot.process()
    elif (action == "MERGEGEM"):
        mergeGem = MergeGem(opts.inFile, opts.out)
        mergeGem.mergeGem()
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
        df = pd.read_csv(self.geneExpFile, sep="\t", quoting=csv.QUOTE_NONE, comment="#")
        if "MIDCounts" in df.columns:
            df.rename(columns={"MIDCounts": "UMICount"}, inplace=True)
        elif 'values' in df.columns:
            df.rename(columns={"values": "UMICount"}, inplace=True)
        elif 'MIDCount' in df.columns:
            df.rename(columns={'MIDCount': 'UMICount'}, inplace=True)
        df['x'] = (df['x']/self.binSize).astype(np.uint32)*self.binSize
        df['y'] = (df['y']/self.binSize).astype(np.uint32)*self.binSize
        df['cell'] = df['x'].astype(str) + "-" + df['y'].astype(str)
        bindf = df['UMICount'].groupby([df['cell'], df['geneID']]).sum()
        cells = set(str(x[0]) for x in bindf.index)
        genes = set(str(x[1]) for x in bindf.index)
        cellsdic = dict(zip(cells, range(0, len(cells))))
        genesdic = dict(zip(genes, range(0, len(genes))))
        rows = [cellsdic[x[0]] for x in bindf.index]
        cols = [genesdic[x[1]] for x in bindf.index]
        expMtx = sparse.csr_matrix((bindf.values, (rows, cols)))

        obs = pd.DataFrame(index = cells)
        var = pd.DataFrame(index = genes)
        adata = anndata.AnnData(X = expMtx, obs = obs, var = var)
        logging.info("ammdata generate finished...")
        return adata, genesdic

    def process(self):
        adata, genesdic = self.bin_stat()
        resultFile = os.path.join(self.outdir, "{0}x{0}_adata.h5ad".format(self.binSize))
        adata.write(resultFile)

class CellCluster():
    def __init__(self, geneExpFile, outFile, binSize):
        self.geneExpFile = geneExpFile
        self.outFile = outFile
        self.binSize = binSize

    def scanpyCluster(self):
        import scanpy as sc
        import anndata

        if (self.geneExpFile.endswith(".h5ad")):
            adata = sc.read_h5ad(self.geneExpFile)
        else:
            df = pd.read_csv(self.geneExpFile, sep="\t", quoting=csv.QUOTE_NONE, comment="#")
            if "MIDCounts" in df.columns:
                df.rename(columns={"MIDCounts": "UMICount"}, inplace=True)
            elif 'values' in df.columns:
                df.rename(columns={"values": "UMICount"}, inplace=True)
            elif 'MIDCount' in df.columns:
                df.rename(columns={"MIDCount": "UMICount"}, inplace=True)
            if 'label' not in df.columns:
                df['x'] = (df['x']/self.binSize).astype(np.uint32)*self.binSize
                df['y'] = (df['y']/self.binSize).astype(np.uint32)*self.binSize
                df['label'] = df['x'].astype(str) + "-" + df['y'].astype(str)
            else:
                labelFile = os.path.join(os.path.dirname(self.geneExpFile), "merge_GetExp_gene_labeled_stat.txt")
                labeldf = pd.read_csv(labelFile, sep="\t")
                labeldict=dict(zip(labeldf['label'], labeldf['x'].astype(str)+"_"+labeldf['y'].astype(str)))
                df.replace({'label': labeldict}, inplace=True)
            bindf = df['UMICount'].groupby([df['label'], df['geneID']]).sum()
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
            adata.write(self.outFile)
            del(df)
            del(bindf)
            adata = sc.read_h5ad(self.outFile) 

        adata.layers['raw_data'] = adata.X

        sc.pp.filter_cells(adata, min_genes=0)
        sc.pp.filter_genes(adata, min_cells=0)
        adata.var['mt'] = adata.var_names.str.decode('utf-8').str.startswith(('mt-', 'MT-'))
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

        adata.raw = adata
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        if adata.var.shape[0] < 2000:
            return 0
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
        adata = adata[:, adata.var.highly_variable]
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, svd_solver="arpack")
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
        sc.tl.tsne(adata)
        sc.tl.umap(adata)
        sc.tl.leiden(adata)
        adata.obs['louvain'] = adata.obs['leiden']
        sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', use_raw=False, n_genes=300, pts=True, layer='raw_data')
        adata.write(self.outFile)

class DnbMatting():
    def __init__(self, geneExpFile, outdir, binSize=50):
        self.geneExpFile = geneExpFile
        self.outdir = outdir
        self.binSize = binSize
        os.makedirs(self.outdir, exist_ok=True)

    def dnb_matting(self):
        import cv2
        expColumnTypes = {"barcode": object, "geneID": 'category', "MIDCount": np.uint32}
        geneDf = pd.read_csv(self.geneExpFile, header=None, names=["barcode", "geneID", "MIDCount"], sep="\t", dtype=expColumnTypes, quoting=csv.QUOTE_NONE)
        geneDf['x'] = geneDf['barcode'].apply(lambda x: int(x.split("_")[-2]))
        geneDf['y'] = geneDf['barcode'].apply(lambda x: int(x.split("_")[-1]))
        geneDf.drop(['barcode'], axis=1, inplace=True)

        #generate bin image
        tempDf = geneDf[['x', 'y', 'MIDCount']].copy()
        tempDf['x'] = tempDf['x'].map(lambda x: int(x/self.binSize))
        tempDf['y'] = tempDf['y'].map(lambda x: int(x/self.binSize))
        binDf = tempDf['MIDCount'].groupby([tempDf['x'], tempDf['y']]).sum().reset_index()
        x1, x2, y1, y2 = binDf['x'].min(), binDf['x'].max(), binDf['y'].min(), binDf['y'].max()
        binDf['x'] = binDf['x'] - x1
        binDf['y'] = binDf['y'] - y1
        sparseMt = sparse.csr_matrix((binDf['MIDCount'].astype(np.uint8), (binDf['y'], binDf['x'])))
        img = sparseMt.toarray()
        
        #median filter by kenel of size 3
        median_filtered = ndimage.median_filter(img, size=3)
        #normalize
        Imin, Imax = median_filtered.min(), median_filtered.max()
        Omin, Omax = 0, 255
        a = float(Omax-Omin)/(Imax-Imin)
        b = Omin - a*Imin
        normalize = a*median_filtered + b

        #process image
        gradx = cv2.Sobel(normalize, ddepth=-1, dx=1, dy=0, ksize=-1)
        grady = cv2.Sobel(normalize, ddepth=-1, dx=0, dy=1, ksize=-1)
        gradient = cv2.subtract(gradx, grady)
        gradient = cv2.convertScaleAbs(gradient)
        blurred = cv2.blur(gradient, (3, 3))
        (_, thresh) = cv2.threshold(blurred, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
        kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (25, 25))
        closed = cv2.morphologyEx(thresh, cv2.MORPH_CLOSE, kernel)
        closed = cv2.erode(closed, None, iterations=4)
        closed = cv2.dilate(closed, None, iterations=4)
        (cnts, _) = cv2.findContours(closed.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

        imageFile = os.path.join(self.outdir, "{0}x{0}_image.png".format(self.binSize))
        if(len(cnts)<1):
            filterGene = geneDf
            cv2.imwrite(imageFile, img)
        else:
            c = sorted(cnts, key=cv2.contourArea, reverse=True)[0]
            rect = cv2.minAreaRect(c)
            box = np.int0(cv2.boxPoints(rect))
            Xs = box[...,0]
            Ys = box[...,1]
            bx1, bx2, by1, by2 = (min(Xs)+x1-0.5)*self.binSize, (max(Xs)+x1+0.5)*self.binSize, (min(Ys)+y1-0.5)*self.binSize, (max(Ys)+y1+0.5)*self.binSize
            filterGene = geneDf.loc[(geneDf['x']>=bx1)&(geneDf['x']<=bx2)&(geneDf['y']>=by1)&(geneDf['y']<=by2)]
            cv2.drawContours(img, [box], -1, (255, 0, 0), 1)
            cv2.imwrite(imageFile, img)

        filterGenefile = os.path.join(self.outdir, "merge_GetExp_gene.txt")
        filterGene[['geneID', 'x', 'y', 'MIDCount']].to_csv(filterGenefile, index=None, sep="\t")


class Visualization():
    def __init__(self, geneExpFile, outdir, maxBinSize=200, progress=4, img_out=True):
        self.geneExpFile = geneExpFile
        self.outdir = outdir
        self.maxBinSize = maxBinSize
        self.progress = progress
        self.img_out = img_out
        self.geneDf = pd.read_csv(self.geneExpFile, sep="\t", quoting=csv.QUOTE_NONE)
        if "MIDCounts" in self.geneDf.columns:
            self.geneDf.rename(columns={"MIDCounts": "values"}, inplace=True)
        elif 'UMICount' in self.geneDf.columns:
            self.geneDf.rename(columns = {"UMICount": "values"}, inplace=True)
        elif 'MIDCount' in self.geneDf.columns:
            self.geneDf.rename(columns={"MIDCount": "values"}, inplace=True)
        dtypes = {"x": np.uint32, "y": np.uint32, "geneID": object, "values": np.uint32}
        self.geneDf.astype(dtypes, copy=False)
        self.geneDf['geneID'] = self.geneDf['geneID'].apply(lambda x: str(x).replace("/", "_"))
        os.makedirs(self.outdir, exist_ok=True)
        self.dnbOutdir = os.path.join(outdir, "dnb_merge")
        os.makedirs(self.dnbOutdir, exist_ok=True)

    def write_h5_total(self):
        total_outf = os.path.join(self.outdir, 'stereomics_total.h5')
        hdf5_fh = h5py.File(total_outf, "w")
        hdf5_dnb_group = hdf5_fh.create_group('dnb_merge')
        hdf5_gene_group = hdf5_fh.create_group('gene_merge')

        # save gene table
        gene_table_df = self.geneDf['values'].groupby(self.geneDf['geneID']).sum().reset_index()
        gene_table_df = gene_table_df.sort_values(by=['values'], ascending=False)
        gene_table_group = hdf5_gene_group.create_group('gene_table')
        gene_list = gene_table_df['geneID'].astype('S')
        gene_table_group['Gene'] = gene_list
        gene_table_group.create_dataset('MIDCounts', data=gene_table_df['values'], dtype='int32')

        # save dnb range
        dnb_range_dict = {'min_x': int(self.geneDf['x'].min()), 'max_x':int(self.geneDf['x'].max()), 'min_y':int(self.geneDf['y'].min()), 'max_y':int(self.geneDf['y'].max())}
        dt = h5py.special_dtype(vlen=str)
        dnb_range_arr = hdf5_dnb_group.create_dataset('dnb_range', (1,), dtype=dt)
        dnb_range_arr[0] = json.dumps(dnb_range_dict)

        hdf5_fh.close()

    def write_bin_h5(self, geneBinDf, bin_size, img_out):
        bin_file_name = 'stereomics_' + str(bin_size) +  '.h5'
        bin_outf = os.path.join(self.outdir, bin_file_name)
        hdf5_fh_bin = h5py.File(bin_outf, "w")
        hdf5_dnb_group_bin = hdf5_fh_bin.create_group('dnb_merge')
        hdf5_gene_group_bin = hdf5_fh_bin.create_group('gene_merge')

        ##gene
        merge_gene_dff = self.merge_gene_v2(geneBinDf, bin_size, bin_size)
        h5_gene_bin_group = hdf5_gene_group_bin.create_group(f'bin{bin_size}')
        for gene, value in merge_gene_dff.groupby(merge_gene_dff.geneID):
            h5_gene_bin_group.create_dataset(gene, data=value[['x','y','values']], dtype='int32')

        #dnb
        merge_dnb_dff = self.merge_dnb_v2(merge_gene_dff, bin_size, bin_size)
        if bin_size==200 and img_out:
            if not os.path.exists(self.dnbOutdir):
                os.makedirs(self.dnbOutdir)
            self.getFig(merge_dnb_dff, os.path.join(self.dnbOutdir, 'bin200.png'), scale=20, dpi=72)
        hdf5_dnb_group_bin.create_dataset(f'bin{bin_size}', data=merge_dnb_dff, dtype='int32')

        hdf5_fh_bin.close()


    def merge_gene_v2(self, gene_df, dot_x, dot_y):
        #gene merge
        gene_dff = gene_df.copy()
        if (dot_x > 1 or dot_y > 1):
            gene_dff['x'] = (gene_dff['x']/dot_x).astype('int')*dot_x
            gene_dff['y'] = (gene_dff['y']/dot_y).astype('int')*dot_y
        gene_dff = gene_dff['values'].groupby([gene_dff['x'], gene_dff['y'], gene_dff['geneID']]).sum().reset_index()
        return gene_dff

    def merge_dnb_v2(self, dnb_df, dot_x, dot_y):
        dnb_dff = dnb_df.copy()
        #dnb_dff['x'] = (dnb_dff['x']/dot_x).astype('int')*dot_x
        #dnb_dff['y'] = (dnb_dff['y']/dot_y).astype('int')*dot_y
        dnb_dff_tmp = dnb_dff['values'].groupby([dnb_dff['x'], dnb_dff['y']]).sum().reset_index()
        gene_count = dnb_dff['geneID'].groupby([dnb_dff['x'], dnb_dff['y']]).nunique().reset_index()
        dnb_dff_tmp['gene_counts'] = gene_count['geneID']
        return dnb_dff_tmp

    def getFig(self, data, outfile, scale=1, dpi=72):
        try:
            cmap = mpl.colors.ListedColormap(['#0C3383', '#0A88BA', '#F2D338', '#F28F38', '#D91E1E'])
            x_range=max(data['x']) - min(data['x'])
            y_range=max(data['y']) - min(data['y'])
            x_num = len(data['x'].drop_duplicates())

            plt.figure(figsize=(1*scale,y_range/x_range*scale), facecolor='#262B3D', edgecolor='black') ## 设置图像大小 inch
            ##去掉图像旁边的空白区
            plt.gca().xaxis.set_major_locator(plt.NullLocator())
            plt.gca().yaxis.set_major_locator(plt.NullLocator())
            plt.gca().xaxis.set_ticks_position('top')
            plt.gca().invert_yaxis()
            plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
            plt.margins(0,0)
            ##添加标题
            #plt.title('Interesting Graph',loc ='center')

            x=data['x']
            y=data['y']
            color=data['values']
            #factor = math.ceil(scale/5)
            #dot_size = math.ceil((dpi * scale)*factor/x_num)
            r = scale*72/(x_range/200)
            dot_size = r**2
            plt.scatter(x, y, c=color, s=dot_size, cmap=cmap)
            plt.axis('off')
            plt.savefig(outfile,facecolor='#262B3D', dpi=dpi, pad_inches = 0)
            return True
        except Exception as e:
            print(e)
            return False

    def copy(self, h5_in, h5_out):
        if 'gene_merge' in h5_in.keys():
            for i in h5_in['gene_merge']:
                if i not in h5_out['gene_merge'].keys():
                    h5_in.copy('gene_merge/' + i, h5_out['gene_merge'])
        if 'dnb_merge' in h5_in.keys():
            for i in h5_in['dnb_merge']:
                if i not in h5_out['dnb_merge'].keys():
                    h5_in.copy('dnb_merge/' + i, h5_out['dnb_merge'])

    def h5_join(self):
        d_names = os.listdir(self.outdir)
        final_outf = os.path.join(self.outdir, 'stereomics.h5')
        h5_out = h5py.File(final_outf, "w")
        h5_out.create_group('gene_merge')
        h5_out.create_group('dnb_merge')
        for h5name in d_names:
            if h5name.endswith('h5') and h5name != 'stereomics.h5':
                full_file_name = os.path.join(self.outdir, h5name)
                h5_in = h5py.File(full_file_name, mode = 'r')
                self.copy(h5_in, h5_out)
                h5_in.close()
                os.remove(full_file_name)

        h5_out.close()

    def process(self):
        binSizeList = [1,2,5,10,15,20,50,80,100,150,200]
        binSizeList = filter(lambda x: x<=self.maxBinSize, binSizeList)
        self.write_h5_total()
       
        if (self.progress == 1):
            for binSize in binSizeList:
                self.write_bin_h5(self.geneDf, binSize, self.img_out)
        else:
            pool = Pool(self.progress)
            for binSize in binSizeList:
                pool.apply_async(self.write_bin_h5, (self.geneDf, binSize, self.img_out,))
            pool.close()
            pool.join()
        self.h5_join()

class ConvertBinData():
    """
    input: The lasso bin gene expression matrix; The complete gene expression matrix
    return: Binsize=1 gene expression matrix.
    """

    def __init__(self, partfile, genefile, outfile, binSize):
        self.typeColumn = {"geneID": 'str', "x": np.uint32, \
                            "y": np.uint32, "values": np.uint32, 'MIDCount':np.uint32, \
                            "MIDCounts":np.uint32, "UMICount": np.uint32}
        self.partfile = partfile
        self.genefile = genefile
        self.outfile = outfile
        self.binSize = binSize

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
        bindf['x'] = df['x'].apply(lambda x: x-self.Xmin if x>self.Xmin else 0)
        bindf['y'] = df['y'].apply(lambda x: x-self.Ymin if x>self.Ymin else 0)
        #bindf['x'] = df['x'] - self.Xmin
        #bindf['y'] = df['y'] - self.Ymin
        bindf['values'] = [255] * len(df)
        
        sparseMt = sparse.csr_matrix((bindf['values'].astype(np.uint8), (bindf['y'], bindf['x'])))
        img = sparseMt.toarray()
        return img

    def ConvertData(self):
        import cv2
        ### Initial data
        if self.binSize > 50:
            raise ValueError("Binsize could not larger than 50.")

        genedf = pd.read_csv(self.genefile, sep='\t', dtype=self.typeColumn, comment="#")
        partdf = pd.read_csv(self.partfile, sep='\t', dtype=self.typeColumn, comment="#")
        
        if len(genedf) < len(partdf):
            raise Warning("Inputs are not correct.")
        self.Xmin = genedf['x'].min()
        self.Ymin = genedf['y'].min()
        partdf = partdf.loc[(partdf['x']>=self.Xmin)&(partdf['y']>=self.Ymin)]
        ### Create Mask
        part_img = self.__CreateImg(partdf)
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2*self.binSize, 2*self.binSize))
        mask = cv2.morphologyEx(part_img, cv2.MORPH_CLOSE, kernel)

        ### Dump results
        mergedf = self.__Dumpresult(mask, genedf)

        mergedf.to_csv(self.outfile, sep="\t", index=None)

class SaturationPlot():
    def __init__(self, saturationStat, outdir):
        self.saturationStat = saturationStat
        self.outdir = outdir
    def process(self):
        import matplotlib.pyplot as plt 
        os.makedirs(self.outdir, exist_ok=True)
        sadf = pd.read_csv(self.saturationStat, sep=" ", quoting=csv.QUOTE_NONE)
        fig = plt.figure(figsize=(12,5),dpi=100)
        ax = plt.subplot(1, 2, 1)
        ax.plot(sadf['bar_x'], sadf['bar_y1'])
        ax.set_xlabel("Total reads number of sampling")
        ax.set_ylabel("Sequencing Saturation")
        ax.grid()
        ax = plt.subplot(1, 2, 2)
        ax.plot(sadf['bar_x'], sadf['bar_y2'])
        ax.set_xlabel("Total reads number of sampling")
        ax.set_ylabel("Median Genes per barcode")
        ax.grid()
        figFilePath = os.path.join(self.outdir, "plot_1x1_saturation.png")
        plt.savefig(figFilePath, format="png")
        fig=plt.figure(figsize=(10,4),dpi=100)
        plt.clf()
        ax = plt.subplot(1,2,1)
        ax.plot(sadf['bin_x'], sadf['bar_y1'])
        ax.set_xlabel("Total reads number of sampling")
        ax.set_ylabel("Sequencing Saturation")
        ax.grid()
        ax = plt.subplot(1, 2, 2)
        ax.plot(sadf['bin_x'], sadf['bin_y2'])
        ax.set_xlabel("Total reads number of sampling")
        ax.set_ylabel("Median Genes per bin")
        ax.grid()
        figFilePath = os.path.join(self.outdir, "plot_150x150_saturation.png")
        plt.savefig(figFilePath, format="png")

class MergeGem():
    def __init__(self, infile, outfile):
        self.infiles = infile.strip().split(",")
        self.outfile = outfile
        self.typeColumn = {"geneID": 'str', "x": np.uint32, \
                            "y": np.uint32, "values": np.uint32, 'MIDCount':np.uint32, \
                            "MIDCounts":np.uint32, "UMICount": np.uint32}
        self.genedf = pd.DataFrame(columns= ['geneID', 'x', 'y', 'MIDCount'])
    def mergeGem(self):
        for gemfile in self.infiles:
            df = pd.read_csv(gemfile, sep="\t", dtype=self.typeColumn, quoting=csv.QUOTE_NONE, comment="#")
            if "MIDCounts" in df.columns:
                df.rename(columns={"MIDCounts": "MIDCount"}, inplace=True)
            elif 'values' in df.columns:
                df.rename(columns={"values": "MIDCount"}, inplace=True)
            elif 'UMICount' in df.columns:
                df.rename(columns={'UMICount': 'MIDCount'}, inplace=True)
            if self.genedf.empty:
                self.genedf = df
            else:
                self.genedf.append(df, ignore_index=True)
        self.genedf.to_csv(self.outfile, sep="\t")

if __name__ == "__main__":
    main()
