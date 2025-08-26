import numpy as np
import scipy as sp
from scipy import stats
import matplotlib.pyplot as plt
import os
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score
import math

def changeIDS(direct):
    data = pd.read_csv(direct, sep="\s+",index_col=False)
    data['FID'] = data['FID'].str.split('_').str[0]
    data['IID'] = data['IID'].str.split('_').str[0]
    data.to_csv(direct,sep="\t",index=False)
# Modify the bim file for the test data.
bimfile = pd.read_csv("test.bim",header=None,sep="\s+")
bimfile[0] = 21
bimfile[2] = list(range(1,len(bimfile)+1))
bimfile.to_csv("test.bim",header=False, index=False, sep="\t")
# Modify the GWAS file, and use the correct chromosome number for each SNP.
data = pd.read_csv("Data.txt",sep="\s+")
data['CHR']=21
data.to_csv("Data.txt.gz", index=False, sep="\t",compression='gzip')
changeIDS("testYRI.covariate")
changeIDS("testYRI.pheno")
changeIDS("trainYRI.covariate")
changeIDS("trainYRI.pheno")
