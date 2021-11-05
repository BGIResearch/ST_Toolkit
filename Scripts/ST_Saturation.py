# @Author: LiuXing liuxing2@genomics.cn 
# @Date: 2021-11-03 17:14:52 
# @Last Modified by:   LiuXing 
# @Last Modified time: 2021-11-03 17:14:52 

import os, sys
import pandas as pd
import numpy as np
from scipy import sparse
import csv
import gzip
import saturation
from optparse import OptionParser
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm

def main():
    """
    %prog [options]
    This softeware can be used to do saturation statistic
    """

    parser = OptionParser(main.__doc__)
    parser.add_option("-i", "--inFile", action="store", type="str", dest="inFile", help="input gene expression matrix contains reads count information")
    parser.add_option("-o", "--out", action="store", type="str", dest="out", help="output directory")
    parser.add_option("--tissue", action="store", type="str", dest="tissue", help="gene expression matrix only contians data under the tissue region")

    opts, args = parser.parse_args()
    if (opts.inFile == None or opts.tissue == None):
        sys.exit(not parser.print_help())

    saturationFile = os.path.join(opts.out, "sequence_saturation.txt")
    os.makedirs(opts.out, exist_ok=True)
    getSaturationTable(opts.inFile, opts.tissue, saturationFile)
    getSaturationFig(saturationFile, opts.out)

def getSaturationTable(rawGem, tissueGem, outFile):
    offsetX, offsetY = _getOffset(tissueGem)
    expColumnTypes = {"x": np.uint32, "y": np.uint32, "geneID": 'category', "MIDCount": np.uint32, "UMI": np.uint64, "reads": np.uint16}
    tissuedf = pd.read_csv(tissueGem, sep="\t", comment="#", dtype=expColumnTypes, usecols=['x', 'y'], quoting=csv.QUOTE_NONE)
    tissuedf['x'] = tissuedf['x'] + offsetX
    tissuedf['y'] = tissuedf['y'] + offsetY
    tissuedf.drop_duplicates(inplace=True)
    coor = set([x[1] << 32 + x[2] for x in tissuedf.itertuples()])
    saturation.saturation(rawGem, outFile, coor)

def getSaturationFig(saturationFile, outdir, binSize=200, readsScale=1000000):
    os.makedirs(outdir, exist_ok=True)
    sadf = pd.read_csv(saturationFile, sep=" ", quoting=csv.QUOTE_NONE)
    fig = plt.figure(figsize=(12,5),dpi=100)
    ax = plt.subplot(1, 2, 1)
    requiredSa = 0.8
    ax.plot(sadf['bar_x']/readsScale, sadf['bar_y1'])
    ax.set_xlabel("Total reads number of sampling (M)")
    ax.set_ylabel("Sequencing Saturation")
    plt.axhline(requiredSa, color='green', lw=1, alpha=0.7)
    ax.grid()
    ax = plt.subplot(1, 2, 2)
    ax.plot(sadf['bar_x']/readsScale, sadf['bar_y2'])
    ax.set_xlabel("Total reads number of sampling (M)")
    ax.set_ylabel("Median Genes per barcode")
    ax.grid()
    figFilePath = os.path.join(outdir, "plot_1x1_saturation.png")
    plt.tight_layout()
    plt.savefig(figFilePath, format="png", bbox_inches='tight')
    plt.clf()
    
    #plot saturation figure of bin200
    fig=plt.figure(figsize=(16,5),dpi=100)
    ax = plt.subplot(1,3,1)
    ax.plot(sadf['bin_x']/readsScale, sadf['bar_y1'])
    ax.set_xlabel("Total reads number of sampling (M)")
    ax.set_ylabel("Sequencing Saturation")
    plt.axhline(requiredSa, color='green', lw=1, alpha=0.7)
    ax.grid()
    ax = plt.subplot(1, 3, 2)
    ax.plot(sadf['bin_x']/readsScale, sadf['bin_y2'])
    ax.set_xlabel("Total reads number of sampling (M)")
    ax.set_ylabel("Median Genes per bin")
    ax.grid()

    xData = sadf['bin_x']
    yData = sadf['bin_umi']
    Rsquared, fittedParameters = _cur_fit(xData, yData)
    threshold = 5000
    maxX = _reFunc(threshold, *fittedParameters)
    ax = plt.subplot(1, 3, 3)
    ax.plot(xData/readsScale, yData, marker='o')
    xModel = np.linspace(min(xData), max(xData))
    yModel = _func(xModel, *fittedParameters)
    xModel1 = np.linspace(max(xData), maxX)
    yModel1 = _func(xModel1, *fittedParameters)
    ax.plot(xModel/readsScale, yModel, color="red", lw=1)
    bstr=str(round(fittedParameters[1], 2)) if fittedParameters[1]<0 else "+{0}".format(round(fittedParameters[1], 2))
    cstr=str(round(fittedParameters[2], 2)) if fittedParameters[2]<0 else "+{0}".format(round(fittedParameters[2], 2))
    rstr="R\u00b2={:0.3f}".format(Rsquared)
    labelstr='y={0}*log(x{1}){2}\n{3}'.format(round(fittedParameters[0],2), bstr, cstr, rstr)
    ax.plot(xModel1/readsScale, yModel1, color="red", linestyle='dashed', label=labelstr, lw=0.5)
    ax.set_xlabel("Total reads number of sampling (M)")
    ax.set_ylabel("unique reads number per bin")
    ax.grid()
    ax.legend()
    plt.axhline(threshold, color='green', lw=1, alpha=0.7)
    plt.axvline(maxX/readsScale, color='green', lw=1, alpha=0.7)
    plt.text(maxX/readsScale,threshold,(int(maxX/readsScale),threshold),color='b')
    plt.text(max(xData)/readsScale,max(yData),(int(max(xData)/readsScale),max(yData)),color='b')
    
    figFilePath = os.path.join(outdir, "plot_{0}x{0}_saturation.png".format(binSize))
    plt.tight_layout()
    plt.savefig(figFilePath, format="png", bbox_inches='tight')

def _getOffset(tissueGem):
    offsetX = 0
    offsetY = 0
    reader = gzip.open(tissueGem, 'rt') if tissueGem.endswith(".gz") else open(tissueGem, 'r')
    for line in reader:
        if not line.startswith("#"):
            break
        elif "OffsetX" in line:
            offsetX = int(line.strip().split("=")[1])
        elif "OffsetX" in line:
            offsetY = int(line.strip().split("=")[1])
    return offsetX, offsetY

def _cur_fit(xData, yData):
    initialParameters = np.array([1.0, 1.0, 1.0])
    fittedParameters, pcov = curve_fit(_func, xData, yData, initialParameters)
    modelPredictions = _func(xData, *fittedParameters)
    absError = modelPredictions - yData
    SE = np.square(absError) # squared errors
    MSE = np.mean(SE) # mean squared errors
    RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE
    Rsquared = 1.0 - (np.var(absError) / np.var(yData))
    return Rsquared, fittedParameters

def _func(x, a, b, c):
    result = a*np.log(x+b)+c
    return result

def _reFunc(y, a, b, c):
    result = np.power(np.e, (y-c)/a) - b
    return result

if __name__ == "__main__":
    main()