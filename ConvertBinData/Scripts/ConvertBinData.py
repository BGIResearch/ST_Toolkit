'''
 # @ Author: Xiaoxuan Tang
 # @ Create Time: 2021-01-06 11:26:59
 # @ Modified by: xiaoxuan Tang
 # @ Modified time: 2021-01-06 14:26:05
 # @ Description:
 '''

import sys, os
import cv2
import pandas as pd
import numpy as np
from scipy import sparse

class ConvertBinData():
    """
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

    def ConvertData(self, partfile, genefile, binSize):
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

        return mergedf

