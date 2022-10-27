# #nCVPy: Evaluate robustness of clock genes in population scale data
# #nCVnetFun: a subfunction used to test whether there is a functional clock network in population scale data.
import re
import sys
import random
import pandas as pd
import mantel_test as mt


def nCVnetFun(finputd, fbenchd, hm="human", nperm=1000, seedN=10):
    # check whether finpud and fbenchd are dataframe
    if isinstance(finputd, pd.DataFrame) and isinstance(fbenchd, pd.DataFrame):
        if (hm.strip()).lower() == "human":       # strip remove the whitespace at beginning and the ending
            fbenchd.iloc[:, 0] = [sig.upper() for sig in list(fbenchd.iloc[:, 0])]
            fbenchd.columns = [sig.upper() for sig in fbenchd.columns]
        # replace the column name
        inputdFirstColumnName = finputd.columns[0]
        finputd = finputd.rename(columns={inputdFirstColumnName: "geneSym"})
        benchdFirstColumnName = fbenchd.columns[0]
        fbenchd = fbenchd.rename(columns={benchdFirstColumnName: "geneSym"})
        # check whether there is duplicate gene symbols
        inputdFirstColumn = list(finputd.loc[:, "geneSym"])
        benchdFirstColumn = list(fbenchd.loc[:, "geneSym"])
        if (len(inputdFirstColumn) == len(set(inputdFirstColumn))) and \
                (len(benchdFirstColumn) == len(set(benchdFirstColumn))):
            fbenchd.index = fbenchd.loc[:, "geneSym"]
            finputd.index = finputd.loc[:, "geneSym"]
            # get the overlapped genes
            bothID = set(fbenchd.index).intersection(finputd.index)
            if len(bothID) >= 3:
                # the reference correlation matrix
                corR = fbenchd.loc[bothID, bothID]
                # the correlation matrix of finputd
                tinputd = finputd.loc[bothID, ]
                tinputd = tinputd.iloc[:, 1:].transpose()
                corD = tinputd.corr(method='spearman')
                corD.index = bothID
                corD.columns = bothID
                # calculate the similarity between two matrices
                simD = mt.mant_test(corD, corR, nperm=1)
                # start the permutation step
                random.seed(seedN)
                indexM = pd.DataFrame([])
                for np in list(range(nperm)):
                    addIndex = random.sample(range(finputd.shape[0]), len(bothID))
                    if (np ==0):
                        indexM.insert(loc=np, column=np, value=addIndex)
                    else:
                        tepIndexM = pd.DataFrame({np: addIndex})
                        indexM = pd.concat([indexM, tepIndexM], axis=1)

                def applyFun(rowIndex, rawInputd, benchCorD):
                    permInputd = rawInputd.iloc[list(rowIndex), :]
                    permInputd = permInputd.iloc[:, 1:].transpose()
                    permcorD = permInputd.corr(method='spearman')
                    permZstat = mt.mant_zstat(permcorD, benchCorD)
                    return permZstat
                pstat = indexM.apply(applyFun, rawInputd=finputd, benchCorD=corR, axis="index")
                pstat = list(pstat)
                realz = simD['z.stat']
                greatNum = [num for num in pstat if num > realz]
                pva = len(greatNum)/nperm
                pvaOut = {"zstat": pd.DataFrame({"tag": ["nCVnet"], "geneNumber": [len(bothID)],
                                                 "zstat": [simD['z.stat']], "pvalue": [pva]}),
                          "npermV": pstat, "cmatrix": corD}
                return pvaOut
            else:
                print("Less than 3 genes are overlapped with the bench gene list. Please check the input data.\n")
        else:
            print("Duplicate gene symbols found in 'finputd' or 'fbenchd'. Please remove duplicated gene symbols.\n")
    else:
        print("Please make sure that 'finputd' and 'fbenchd' are data frame.\n")


if __name__ == "__main__":
    import time
    run_start = time.time()
    if len(sys.argv) < 4:
        print("Usage: python nCVnet.py datafile benchfile outfile "
              "human_mouse(optional) nperm(optional) seedN(optional)\n")
        sys.exit(0)
    # check whether the datafile is csv or txt
    if re.search(r'\.csv', sys.argv[1]):
        inputd = pd.read_csv(sys.argv[1])
    elif re.search(r'\.txt', sys.argv[1]):
        inputd = pd.read_table(sys.argv[1])
    else:
        print("the 'datafile' should be a csv or txt file")
        sys.exit(0)
    # check whether the benchfile is csv or txt
    if re.search(r'\.csv', sys.argv[2]):
        benchd = pd.read_csv(sys.argv[2])
    elif re.search(r'\.txt', sys.argv[2]):
        benchd = pd.read_table(sys.argv[2])
    else:
        print("the 'benchfile' should be a csv or txt file")
        sys.exit(0)
    # check other parameters
    ihm = "human"
    inperm = 1000
    iseedN = 10
    if len(sys.argv) >= 5:
        ihm = sys.argv[4]
    if len(sys.argv) >= 6:
        inperm = int(sys.argv[5])
    if len(sys.argv) >= 7:
        iseedN = int(sys.argv[6])
    # run the nCVnet
    netvd = nCVnetFun(inputd, benchd, hm=ihm, nperm=inperm, seedN=iseedN)
    netvd = netvd["zstat"]
    # output the dataframe
    netvd.to_csv(sys.argv[3], index=False)
