import pandas as pd
import os
from sklearn.model_selection import train_test_split
import sys

def saveandformatsamplefile(top,bottom,direc):
### Split test cases/controls and train cases/controls
###make a phenotype file
#bottom['ID_2'] = 'sample_' + bottom['ID_2'].astype(str)
#bottom['ID_1'] = 'sample_' + bottom['ID_1'].astype(str)
    phenotype = pd.DataFrame()
    phenotype["user_id"] = bottom["ID_1"].values
    phenotype["phenotype"] = bottom["pheno"].values
    phenotype.to_csv("phenotype.csv",sep="\t",index=False)
###make covarite file
    cov = pd.DataFrame()
    cov["FID"] = bottom["ID_2"].values
    cov["IID"] = bottom["ID_1"].values
    cov["Sex"] = 1
    cov["cov1"] = bottom["sharedConfounder1_bin1"].values
    cov["cov2"] = bottom["independentConfounder2_cat_norm1"].values
    cov["cov3"] = bottom["independentConfounder2_cat_norm2"].values
    cov["cov4"] = bottom["independentConfounder3_cat_unif1"].values
    sampletop = top.copy()
    samplebottom = bottom.copy()
    samplebottom["pheno"] = samplebottom["pheno"].apply(pd.to_numeric)
    samplebottom.pheno[samplebottom['pheno']<0]=0
    samplebottom.pheno[samplebottom['pheno']>0]=1
    samplebottom["pheno"] = pd.to_numeric(samplebottom["pheno"],downcast='integer')
    sample = pd.concat([sampletop, samplebottom], axis=0)
    sample= sample.astype(str)
    if "test" in direc:
        subsubsection(phenotype,direc,"test")
        sample.to_csv("test_snptest.sample",index=False,sep=" ")
    if "train" in direc:
        subsubsection(phenotype,direc,"train")
        sample.to_csv("train_snptest.sample",index=False,sep=" ")
# Modify the cases/controls information because plink considers 1 as a control and 2 
# as a case, whereas other tools consider 0 as control and 1 as a case.
    sample.pheno[sample['pheno']=='1']='2'
    sample.pheno[sample['pheno']=='0']='1'
    data = sample[["ID_1","ID_2","missing","pheno"]]
    if "test" in direc:
        data.to_csv("test.sample",index=False,sep=" ")
    if "train" in direc:
        data.to_csv("train.sample",index=False,sep=" ")
    samplebottom.pheno[samplebottom['pheno']==1]='2'
    samplebottom.pheno[samplebottom['pheno']==0]='1'

    cov["cov5"] = bottom["sharedConfounder4_norm1"].values
    cov["cov6"] = bottom["independentConfounder4_norm1"].values
    cov.to_csv(direc+"YRI.covariate",index=False,sep="\t")
###PRS phenotype
    phenotype = pd.DataFrame()
    phenotype["FID"] = bottom["ID_2"].values
    phenotype["IID"] = bottom["ID_1"].values
    phenotype["phenotype"] = bottom["pheno"].values
    phenotype.to_csv(direc+"YRI.pheno",sep="\t",index=False)
###NewSample file
    sample = pd.concat([top, bottom], axis=0)
    sample = sample[['ID_1','ID_2', 'missing','pheno']]
    sample.to_csv("YRIs.sample",index=False,sep=" ")
    return phenotype['FID'].values

def splitsample(sample,direc):
# Extract the first row because it does not contain the sample information.
    sampletop = sample.head(1)
    samplebottom = sample.tail(len(sample)-1)
#Modify the sample ID's.
    samplebottom['ID_1'] = samplebottom['ID_1'].astype(str)+str("_") + samplebottom['ID_1'].astype(str)
    samplebottom['ID_2'] = samplebottom['ID_2'].astype(str)+str("_") + samplebottom['ID_2'].astype(str)
    #samplebottom['ID_1'] = samplebottom['ID_1'].astype(str)
    #samplebottom['ID_2'] = samplebottom['ID_2'].astype(str)
    samplebottom["pheno"] = samplebottom["pheno"].apply(pd.to_numeric)
    #PhenotypeSimulator generates continuous phenotype, which we converted to binary phenotype by thresholding on 0.
    samplebottom["pheno"].values[samplebottom["pheno"] < 0] = 0
    samplebottom["pheno"].values[samplebottom["pheno"] > 0] = 1
    samplebottom["pheno"] = pd.to_numeric(samplebottom["pheno"],downcast='integer')
# Spit the samples. The default is 75 percent training and 25 percent test sets.
    x_train, x_test, y_train, y_test = train_test_split(samplebottom, samplebottom["pheno"].values)

    sampletop.iloc[0,9]="B"
    trainsample = saveandformatsamplefile(sampletop, x_train,"train")
    testsample = saveandformatsamplefile(sampletop, x_test,"test")
    return trainsample,testsample

def commit(name):
    sample = pd.read_csv(name+".sample",sep=" ")
    samplebottom = sample.tail(len(sample)-1)
    fam = pd.read_csv(name+".fam",sep="\s+",header=None)
    fam[5] = samplebottom['pheno'].values
    fam.to_csv(name+".fam",header=False,index=False, sep=" ")

# Read the sample files, and ensure path is correct.
originalsamples = pd.read_csv("Ysim_snptest.sample",sep=" ")
# This function splits the samples into training and test sets.
train,test = splitsample(originalsamples,'train')
train,test = splitsample(originalsamples,'test')
