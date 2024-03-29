import sys,os
import random
import numpy as np

BIN_SIZE = 200
SAMPLE_RATIOS = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]

# Read raw data from file, and filtered by coordinates
# raw data format: y x gene umi cnt
def _getData(filename, uniqCoordinates):
    res = []
    uniqBinBarcodes = set()
    with open(filename) as fh_in:
        for line in fh_in:
            line = line.strip().split()
            barcode = (int(line[1]) << 32) + int(line[0])
            if uniqCoordinates is None or barcode in uniqCoordinates:
                res.append(list(map(int, line)))
                uniqBinBarcodes.add((int(line[0])//BIN_SIZE << 32) + int(line[1])//BIN_SIZE)
    return res,uniqBinBarcodes

# Calculate sequencing saturation and medium gene numbers
def _calculate(nestMapData, realReads, isBin=False, umisData=None):
    n_reads = 0
    n_uniq = 0
    n_genes = []
    n_umis = []
    for v1 in nestMapData.values():
        genes = set()
        umis = set()
        for key,v2 in v1.items():
            n_reads += v2
            gene = (key >> 64)
            genes.add(gene)
            umi = key & 0xFFFFFFFFFFFFFFFF
            umis.add(umi)
        n_uniq += len(v1)
        n_genes.append(len(genes))
        n_umis.append(len(umis))
    res = ""
    ratio = 0.0 if n_reads == 0 else 1-(n_uniq*1.0/n_reads)
    medianGene = 0 if len(n_genes) == 0 else int(np.median(n_genes))
    if isBin:
        umis = []
        for d in umisData:
            umis.append(len(umisData[d]))
        meanUmi = 0 if len(umis) == 0 else int(np.mean(umis))
        #meanUmi = 0 if len(n_umis) == 0 else int(np.mean(n_umis))
        res = " {} {:.7} {} {}".format(realReads, ratio, medianGene, meanUmi)
    else:
        if n_reads != 0:
            n_uniq = int(n_uniq*realReads/n_reads)
        res = " {} {:.7} {} {}".format(realReads, ratio, medianGene, n_uniq)
    return res

# Sampling to reduce the amount of data to be calculated
def _sample(rawData, uniqCoordinates, sampleRatio):
    res = []
    sampleNumber = int(len(uniqCoordinates)*sampleRatio)
    # if sampleNumber < 100:
    #     sampleNumber = len(uniqCoordinates)
    sampledCoordinates = set(random.sample(uniqCoordinates, sampleNumber))
    for line in rawData:
        barcode = (line[0]//BIN_SIZE << 32) + line[1]//BIN_SIZE
        if barcode in sampledCoordinates:
            res.extend([line[:-1]] * line[-1])
    return res

def fakeUniqCoordinates(filename):
    res = set()
    with open(filename) as fh_in:
        for line in fh_in:
            line = line.strip().split()
            barcode = (int(line[1]) << 32) + int(line[0])
            res.add(barcode)
    return res

# Main entrance
# inputFile: format like, "y x gene umi count"
# outputFile: result of sequencing saturation and medium gene types
# uniqCoordinates: set of coordinates under tissue area, the element type is int64, means ((x<<32) + y)
# sampleRatio: the percentage of coordinates are taken
def saturation(inputFile, outputFile, uniqCoordinates=None, sampleRatio=0.05):
    data,uniqBinBarcodes = _getData(inputFile, uniqCoordinates)
    if len(data) == 0:
        print("No data leave after filter by coordinates!")
        return

    _saturation(inputFile, outputFile, data, uniqBinBarcodes, uniqCoordinates, sampleRatio)

def _saturation(inputFile, outputFile, data, uniqBinBarcodes, uniqCoordinates=None, sampleRatio=0.05):
    #print("start saturation, unique coordinates numbers", len(uniqCoordinates))
    #print("raw data number {}, unique bin coordinates numbers {}".format(len(data), len(uniqBinBarcodes)))
    totalReadsNum = 0
    for line in data:
        totalReadsNum += line[-1]
    # if uniqCoordinates is not None:
    data = _sample(data, uniqBinBarcodes, sampleRatio)
    #print("data number after sampling", len(data))
    random.shuffle(data)

    pos = 0
    stat_barcode = {}
    stat_bin = {}
    stat_umi = {}
    outStr = "#sample bar_x bar_y1 bar_y2 bar_umi bin_x bin_y1 bin_y2 bin_umi\n"
    for pct in SAMPLE_RATIOS:
        outStr += str(pct)
        sampleEndPos = int(len(data)*pct)
        #print("-----sample:", pct, sampleEndPos)
        realReads = int(totalReadsNum * pct)
        
        while (pos < sampleEndPos):
            (b1, b2, gene, umi) = data[pos]
            key = (gene << 64) + umi
            barcode = (b1 << 32) + b2
            if barcode not in stat_barcode: stat_barcode[barcode] = {}
            if key not in stat_barcode[barcode]: stat_barcode[barcode][key] = 0
            stat_barcode[barcode][key] += 1

            new_barcode = ((b1//BIN_SIZE) << 32) + (b2//BIN_SIZE)
            if new_barcode not in stat_bin: stat_bin[new_barcode] = {}
            if key not in stat_bin[new_barcode]: stat_bin[new_barcode][key] = 0
            stat_bin[new_barcode][key] += 1
            
            if new_barcode not in stat_umi: stat_umi[new_barcode] = set()
            uniq_umi = str(b1)+str(b2)+str(gene)+str(umi)
            stat_umi[new_barcode].add(uniq_umi)
            
            pos += 1

        outStr += _calculate(stat_barcode, realReads, False)
        outStr += _calculate(stat_bin, realReads, True, stat_umi)
        outStr += "\n"

    with open(outputFile, 'w') as fh_out:
        fh_out.write(outStr)
        
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("prog <inputFile> <outputFile>")
        sys.exit(1)

    inputFile = sys.argv[1]
    outputFile = sys.argv[2]
    uniqCoordinates = None
    ## uniqCoordinates = fakeUniqCoordinates(inputFile)
    sampleRatio = 0.05
    saturation(inputFile, outputFile, uniqCoordinates, sampleRatio)
