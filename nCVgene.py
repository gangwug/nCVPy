# ##nCVPy: Evaluate robustness of clock genes in population scale data
# #nCVgeneFun: a subfunction used to calculate nCV values of targeted clock genes.
# #checked the output results with R output results for the example data,
# #the difference is acceptable （< 1e-12),
# #this minor difference is also seen in previous comparison(MetaCycle vs ARSER.py)
# #nCVgene.py is good to use
import re
import sys
import pandas as pd
from statistics import mean, stdev


# #the function of read target gene list
def readTargetGeneList(filename):
    # read seed genes from file, line started with # will be skipped, one gene per line
    tarGeneList = []
    with open(filename) as fin:
        for line in fin:
            if line[0] == "#":
                continue
            elif len(line.strip()) == 0:
                continue
            else:
                tarGeneList.append(line.strip())
    return tarGeneList


# #the function of calculating nCV
def nCVgeneFun(inputd, cgenes="a"):
    cinputd = pd.DataFrame(inputd)
    first_column_name = cinputd.columns[0]
    cinputd = cinputd.rename(columns={first_column_name: "geneSym"})
    first_column = list(cinputd.loc[:, "geneSym"])
    # check whether there is duplicate gene symbols
    if len(first_column) == len(set(first_column)):
        # remove the first column
        rcinputd = cinputd.iloc[:, 1:]
        # calculate the mean value
        rcmeanv = rcinputd.apply(mean, axis=1)
        rcsd = rcinputd.apply(stdev, axis=1)
        rccv = rcsd/rcmeanv
        rcncv = rccv/mean(rccv)
        outd = pd.DataFrame({"geneSym": first_column, "nCV": rcncv, "meanv": rcmeanv})
        # check the cgenes
        cgenes = cgenes.split(";")
        if len(cgenes) > 1:
            # check the intersection
            share_name = set(first_column).intersection(cgenes)
            if len(share_name) > 0:
                in_indice = [cgenes.index(sname) for sname in share_name]
                in_indice.sort()
                cgenes = [cgenes[cindex] for cindex in in_indice]
                out_indice = [first_column.index(cgene) for cgene in cgenes]
                outd = outd.iloc[out_indice, :]
            else:
                outd = pd.DataFrame({"NoOutputA": ["nogene"], "NoOutputB": ["nocv"]})
                print("There is no overlap between 'cgenes' and input data.\n")
    else:
        outd = pd.DataFrame({"NoOutputA": ["nogene"], "NoOutputB": ["nocv"]})
        print("Duplicate gene ids are detected. " +
              "Please merge rows with duplicate gene ids.\n")
    return outd


if __name__ == "__main__":
    import time
    run_start = time.time()
    if len(sys.argv) < 3:
        print("Usage: python nCVgene.py infile outfile targetGeneListFile(optional)")
        sys.exit(0)
    if re.search(r'\.csv', sys.argv[1]):
        inputd = pd.read_csv(sys.argv[1])
    elif re.search(r'\.txt', sys.argv[1]):
        inputd = pd.read_table(sys.argv[1])
    else:
        print("the 'infile' should be a csv or txt file")
        sys.exit(0)
    # check whether there is a target gene list
    icgenes = "None"
    if len(sys.argv) > 3:
        oricgenes = readTargetGeneList(filename=sys.argv[3])
        icgenes = ";".join(oricgenes)
        # print("\,".join(cgenes)+"\n")
    # calculate the nCV
    ncvd = nCVgeneFun(inputd, cgenes=icgenes)
    # output the data
    ncvd.to_csv(sys.argv[2], index=False)
