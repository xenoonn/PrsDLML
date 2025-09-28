import pandas as pd

snpteststats = pd.read_csv("train.sum",sep=" ",low_memory=False,skiprows = 10)
snpteststats = snpteststats.head(len(snpteststats)-1)
plinkstats = pd.read_csv("train.assoc.fisher",sep="\s+",low_memory=False)
gwasstats = pd.DataFrame()
# Change the column names as required by lassosum and plink.
# Ensure chromosome number is correct. We changed it in the later step.
gwasstats['CHR'] = ['0']*len(plinkstats)
gwasstats['BP'] = plinkstats['BP'].values
gwasstats['SNP'] = plinkstats['SNP'].values
gwasstats['A1'] = snpteststats['alleleA'].values
gwasstats['A2'] = snpteststats['alleleB'].values
gwasstats['N'] = snpteststats['all_total'].values
gwasstats['SE'] = snpteststats['frequentist_add_se_1'].values
gwasstats['P'] = plinkstats['P'].values
gwasstats['OR'] = plinkstats['OR'].values
# In simulated data, the imputation score is 1 as we have not used any imputation technique.
gwasstats['INFO'] = 1
gwasstats['MAF'] = snpteststats['all_maf'].values
#gwasstats.to_csv(traindirec+os.sep+"Data.txt.gz",index=False, sep="\t",compression='gzip')
gwasstats.to_csv("Data.txt",index=False, sep="\t")
