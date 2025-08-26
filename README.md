# Simulation of genotype-phenotype for different populations

# Introduction

Polygenic risk scores evaluate the influence of your genes on disease development. To gain a better understanding of your likelihood of developing a specific disease, it is important to integrate polygenic risk scores with other factors that also impact disease risk, rather than relying solely on either one independently. A big amount of features should be taken into account when analyzing PRS, but in our study we pay full attention to populations -  how populations affect polygenic risk for different phenotypes. 

Our goal is to create deep learning algorithms, which can easily and accurately predict polygenic risk for an individual based on his genotype. However, firstly we need to prove our hypothesis - PRS varies on population which an individual has. For this purpose, we simulate data for various populations with possible tools such as *PhenotypeSimulator*, *SimulatePhenotypes* from *HapGen*, *GEPSI*.

# Data Simulation

## Genotype generation with HapGen2

***Hapgen*** is a software tool used for simulating genetic data based on haplotype data. It is commonly used in genetic research to generate synthetic datasets that mimic real-world genetic variation. The process of data simulation from Hapgen involves several steps, which are as follows:

1. Input Data: The first step is to provide the necessary input data to Hapgen. This includes the reference haplotype data, which represents the known genetic variations in a population. It can be obtained from various sources such as the 1000 Genomes Project or other genetic databases. Additionally, Hapgen requires a recombination map, which describes the frequency and location of genetic recombination events.
2. Parameter Specification: Hapgen allows users to specify various parameters to control the simulation process. These parameters include the sample size, which determines the number of individuals in the simulated dataset, and the disease model, which defines the relationship between genetic variants and a specific trait or disease.
3. Haplotype Phasing: Hapgen employs a statistical algorithm to phase the input haplotypes. Phasing involves determining the parental origin of each allele at a given genomic position. This step is crucial for accurately simulating genetic data, as it ensures that the generated genotypes reflect the underlying haplotype structure.
4. Recombination: Hapgen utilizes the provided recombination map to introduce genetic recombination events into the simulated data. Recombination occurs during meiosis when genetic material from parents is shuffled, leading to the creation of new combinations of alleles. By incorporating recombination, Hapgen generates realistic patterns of genetic variation observed in real populations.
5. Genotype Generation: Once haplotype phasing and recombination have been performed, Hapgen generates individual genotypes for each simulated sample. Genotypes represent the specific combination of alleles at each genomic position for an individual. These genotypes can be in various formats, such as VCF (Variant Call Format) or PLINK files, which are commonly used in genetic analysis.
6. Quality Control: After the genotypes are generated, it is essential to perform quality control checks to ensure the integrity of the simulated data. This may involve filtering out low-quality genotypes, removing individuals with excessive missing data, or applying other quality control measures commonly used in genetic studies.
7. Output: Finally, Hapgen produces the simulated genetic dataset, which can be used for various downstream analyses such as association studies, heritability estimation, or testing statistical methods. The output includes the simulated genotypes and any additional information specified during the parameter specification step.

Overall, the process of data simulation from Hapgen involves inputting reference haplotypes and a recombination map, specifying simulation parameters, performing haplotype phasing and recombination, generating individual genotypes, conducting quality control checks, and obtaining the final simulated dataset.

**Here we describe the process of data simulation:**

Download HapGen2 tool:

```python
! wget https://mathgen.stats.ox.ac.uk/genetics_software/hapgen/download/builds/x86_64/v2.2.0/hapgen2_x86_64.tar.gz
! tar -zxvf hapgen2_x86_64.tar.gz
```

First of all, we need to download all necessary data for  simulation. In our study we are using 1000 genomes references panels. It is available here: [https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html](https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html)

```python
! wget https://mathgen.stats.ox.ac.uk/impute/ALL_1000G_phase1integrated_v3_impute.tgz
! tar -xvzf ALL_1000G_phase1integrated_v3_impute.tgz
```

After downloading all the necessary files, we unzip them and get a big amount of files.


The files are a gzipped tar archive and contains 4 kinds of files :

1. ***.hap.gz** (one file per autosome) -- Phased haplotype file in IMPUTE -h format.

2. ***.legend.gz** (one file per autosome) -- Legend file in IMPUTE -l format (compressed by gzip software; variant positions in NCBI b37 coordinates). These files have the following columns:

column 1 (id) - variant ID

column  2 (position) - base pair position

column 3 (a0) - allele labeled '0' in .hap file. This is the REF (or reference) allele.

column 4 (a1) - allele labeled '1' in .hap file. This is the ALT (or alternate) allele.

column 5 (TYPE) - SNP/INDEL/SV denotes type of biallelic variant

column 6-10 (AFR, AMR, EAS, EUR, SAS) - ALT allele frequency in continental groups. The mapping of populations to groups is given below.

column 11 (ALL) - ALT allele frequency across all 2,504 samples

variant IDs  Each bi-allelic variant has an ID of the form rsID:ALT or chromosome:position:ALT.

3. **genetic_map** (one file per autosome) -- Genetic map file in IMPUTE -m format (physical positions in NCBI b37 coordinates).

4. **1000GP_Phase3.sample** -- Text file with sample ID, population and continental group for the individuals in the haplotype files. There are two haplotypes per sample, and the column order of haplotypes matches the row order of sample IDs. The population IDs are defined below, grouped by continental group.

Since every file is heavy weighted, we **settled** on chromosome 21. The files we need are *ALL_1000G_phase1integrated_v3_chr21_impute.hap.gz,*

*ALL_1000G_phase1integrated_v3_chr21_impute.legend.gz,* 

*ALL_1000G_phase1integrated_v3.sample,* 

*genetic_map_chr21_combined_b37.txt*

It is worth mentioning, firstly we decided to split the dataset into three populations - European, African and Asian. As we simulate genotypes for different populations, we need to check coordinates of each population in the sample file. In the file we can see a population and a group. Since we have very limited resources, our simulation is based on a group, not on a population. There are 5 groups in the sample file, we choose only 3 groups - European, Asian and African. To split the 21 chromosome genotype file, we use the coordinates of each group in the sample file and then run this:

For European group:

```bash
cut -f 2-182 -d " " ALL_1000G_phase1integrated_v3_chr21_impute.hap \
> chr21_1.eur_subset.hap
cut -f 359 -d " " ALL_1000G_phase1integrated_v3_chr21_impute.hap \
> chr21_2.eur_subset.hap
cut -f 397-402 -d " " ALL_1000G_phase1integrated_v3_chr21_impute.hap \
> chr21_3.eur_subset.hap
cut -f 405-497 -d " " ALL_1000G_phase1integrated_v3_chr21_impute.hap \
> chr21_4.eur_subset.hap
cut -f 996-1093 -d " " ALL_1000G_phase1integrated_v3_chr21_impute.hap \
> chr21_5.eur_subset.hap

cat chr21_1.eur_subset.hap chr21_2.eur_subset.hap chr21_3.eur_subset.hap chr21_4.eur_subset.hap chr21_5.eur_subset.hap > chr21.eur_subset.hap
```

For African group:

```bash
cut -f 498-995 -d " " ALL_1000G_phase1integrated_v3_chr21_impute.hap \
> chr21.afr_subset.hap
```

For Asian group:

```bash
cut -f 183-226 -d " " ALL_1000G_phase1integrated_v3_chr21_impute.hap \
> chr21_1.asn_subset.hap
cut -f 229-260 -d " " ALL_1000G_phase1integrated_v3_chr21_impute.hap \
> chr21_2.asn_subset.hap
cut -f 265-288 -d " " ALL_1000G_phase1integrated_v3_chr21_impute.hap \
> chr21_3.asn_subset.hap
cut -f 517-613 -d " " ALL_1000G_phase1integrated_v3_chr21_impute.hap \
> chr21_4.asn_subset.hap
cut -f 635-693 -d " " ALL_1000G_phase1integrated_v3_chr21_impute.hap \
> chr21_5.asn_subset.hap
cut -f 702-731 -d " " ALL_1000G_phase1integrated_v3_chr21_impute.hap \
> chr21_6.asn_subset.hap
cat chr21_1.asn_subset.hap chr21_2.asn_subset.hap chr21_3.asn_subset.hap chr21_4.asn_subset.hap chr21_5.asn_subset.hap chr21_6.asn_subset.hap \
> chr21.asn_subset.hap
```

Our next step is running HapGen Generation for each group:

For European group:

```bash
/content/hapgen2 \
-m /content/genetic_map_chr21_combined_b37.txt \
-l /content/ALL_1000G_phase1integrated_v3_chr21_impute.legend \
-h /content/chr21.eur_subset.hap \
-o eur_genotypes_chr21_hapgen \
-n 3000 0 \
-dl -dl 9411618 0 0 0 \
-no_haps_output
```

For African group:

```bash
/content/hapgen2 \
-m /content/genetic_map_chr21_combined_b37.txt \
-l /content/ALL_1000G_phase1integrated_v3_chr21_impute.legend \
-h /content/chr21.afr_subset.hap \
-o afr_genotypes_chr21_hapgen \
-n 3000 0 \
-dl -dl 9411618 0 0 0 \
-no_haps_output
```

For Asian group:

```bash
/content/hapgen2 \
-m /content/genetic_map_chr21_combined_b37.txt \
-l /content/ALL_1000G_phase1integrated_v3_chr21_impute.legend \
-h /content/chr21.asn_subset.hap \
-o asn_genotypes_chr21_hapgen \
-n 3000 0 \
-dl 9411618 0 0 0 \
-no_haps_output
```

Here we describe all flags: 

-m - A file containing the fine-scale recombination rate across the region. 

-l - A legend file for the SNP markers. This file should have 4 columns with one line for each SNP. The columns should contain an ID for each SNP i.e. rs id of the marker, the base pair position of each SNP, base represented by 0 and base represented by 1. The first line of the legend file are column labels (these are not used by the program but the file is required to contain a header line).

-h - File of known haplotypes, with one row per SNP and one column per haplotype. Every haplotype file needs a corresponding legend file (see below), and all alleles must be coded as 0 or 1 -- no other values are allowed.

-o - output file prefix

-n - Sets the number of control and the number of case individuals to simulate. For example **-n 100 200** simulates 100 control and 200 case individuals. The default is to generate 1 control and 1 case individual.

-dl - Sets location, risk allele and relative risks for each disease risk. For each disease SNP, four numbers are required in the following order:

1. physical location of SNP, which must be in the legend file supplied to the -l flag
2. risk allele (0 or 1), the corresponding base can be found in the legend file
3. heterozygote disease risk
4. homozygote disease risk

-no_haps_output - No haplotype data files, ***.haps[.gz]**, will be outputted for the case and control data.

We can see that the only thing that differs is *.hap* file, **the used SNP is the same for every group**. **The disease SNP must not change because, as said before, our main goal is to prove the impact of population in counting PRS.**

After running the commands, as an output we get 5 files, 2 of them are needed for the further work. These files are:

*.controls.gen

*.controls.sample

We generated these files for each population.


## PhenotypeSimulator

***PhenotypeSimulator*** is a software tool used in genetic research to simulate phenotypic data based on genetic and environmental factors. It takes input data such as genetic and environmental variables, and allows users to specify parameters such as heritability and effect sizes. It then applies genetic and environmental effects to generate individual phenotypes, which represent observable traits or characteristics. The simulated dataset can be used for various analyses, such as genetic association studies or heritability estimation. *PhenotypeSimulator* helps researchers mimic real-world phenotypic variation and study complex gene-environment interactions.

So, it’s time to simulate phenotype since we have simulated the genotype.

For phenotype simulation we use the next code:

```bash
library(PhenotypeSimulator)
indir <- "/Users/ico/Desktop/eur"
# specify directory to save data; if it doesn't exist yet, create i.e. needs
# writing permissions
datadir <- '/Users/ico/Desktop/eur/code'
if (!dir.exists(datadir)) dir.create(datadir, recursive=TRUE)
# specify filenames and parameters
genVar <- 0.6
noiseVar <- 1 - genVar
totalSNPeffect <- 0.01
h2s <- totalSNPeffect/genVar
phi <- 0.6
rho <- 0.1
delta <- 0.3
shared <- 0.8
independent <- 1 - shared

#kinshipfile <- paste(indir, "/genotypes_genome_hapgen.controls.grm.rel",
#                     sep="")
genoFilePrefix <- paste(indir, "/genotypes_", sep="")
genoFileSuffix <- "_hapgen.controls.gen"
numberofPeople <- 3000
# simulate phenotype with three phenotype components
phenotype <- runSimulation(N = numberofPeople, P = 1, genoFilePrefix=genoFilePrefix,
                           genoFileSuffix=genoFileSuffix,
                           format = "oxgen", cNrSNP = 5, genVar = genVar, h2s = h2s,
                           phi = 0.6, delta = 0.3, distBetaGenetic = "unif",
                           mBetaGenetic = 0.5, chr=21, 
                           sdBetaGenetic = 1, NrFixedEffects = 4, NrConfounders = c(1, 2,1, 2), 
                           pIndependentConfounders = c(0, 1, 1, 0.5),
                           distConfounders = c("bin", "cat_norm", "cat_unif", "norm"),
                           probConfounders = 0.2, catConfounders = c(3, 4), pcorr = 0.8,
                           verbose = FALSE, seed = 3000)

out <- savePheno(phenotype, directory=datadir,
                 outstring='eur',
                 format=c("csv", "snptest"), verbose=FALSE)
```

As we can observe, many arguments are used in the simulation running.

| **N**  | Number  of samples to simulate |
| --- | --- |
| **P** | Number of phenotypes to simulate |
| **genoFilePrefix**  | Full path/to/chromosome-wise-genotype-file-ending |
| **genoFileSuffix**  | Following chromosome number including .fileformat  |
| **format**  | Name of genotype file format. Options are: "oxgen", "bimbam" or "delim" |
| **cNrSNP**  | Number of causal SNPs; used as genetic variant effects |
| **genVar**  | Proportion of total genetic variance |
| **h2s**  | Proportion of genetic variance of genetic variant effects |
| **phi**  | Proportion of noise variance of observational noise effects; sum of rho, delta and phi has to be equal 1 |
| **delta**  | Proportion of noise variance of non-genetic covariate effects; sum of
rho, delta and phi has to be equal 1 |
| **distBetaGenetic**  | Name of distribution to use to simulate effect sizes of genetic variants; one of "unif" or "norm" |
| **mBetaGenetic**  | Mean/midpoint of normal/uniform distribution for effect sizes of genetic variants |
| **chr**  | Vector of chromosome(s) to chose NrCausalSNPs from |
| **sdBetaGenetic**  | Standard deviation/extension from midpoint of normal/uniform distribution for effect sizes of genetic variants |
| **NrFixedEffects**  | Number of different confounder effects to simulate; allows to simulate fixed effects from different distributions or with different parameters |
| **NrConfounders**  | Vector of number(s) of confounders from a specified distribution to simulate |
| **pIndependentConfounders**  | Vector of proportion(s) of confounders to have a trait-independent effect |
| **distConfounders**  | Vector of name(s) of distribution to use to simulate confounders; one of "unif", "norm", "bin", "cat_norm", "cat_unif" |
| **probConfounders**  | Vector of probability(s) of binomial confounders (0/1); required if distConfounders "bin" |
| **catConfounders**  | Vector of number(s) of confounder categories; required if distConfounders "cat_norm" or "cat_unif" |
| **pcorr**  | Initial strength of correlation between neighbouring traits. Decreases by pcorr^(distance); distance from 0 to P-1 |
| **verbose**  | If TRUE, progress info is printed to standard out |
| **seed**  | Seed to initiate random number generation |

After running all these commands for every group, we get exact files for every group. 

We don’t need all these files, only some of them.

Just in case, our next steps are described here - 

[](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwi1lNGzuZiAAxWAHhAIHV_BDbkQFnoECBkQAQ&url=https://assets.researchsquare.com/files/rs-1298372/v1/893e2612819475eda75a3a4e.pdf&usg=AOvVaw2SPh5gzd8JiEiPt9djJiWa&opi=89978449)

The article:

[An empirical comparison between polygenic risk scores and machine learning for case/control classification](https://www.researchsquare.com/article/rs-1298372/v1)

For further analysis we need to preprocess and structure the resulting data. 

```bash
# Process the final files and save them in VCF file format.
/Users/ico/Desktop/biotools/gtool  -G --g Genotypes_snptest.gen --s \
Ysim_snptest.sample --ped X.ped --map X.map --phenotype pheno
```

```bash
/Users/ico/Desktop/biotools/plink/plink  -file X --recode vcf --out X
```

```python
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
```

```bash
bcftools view -S test_id.txt X.vcf > test.vcf
```

```bash
/Users/ico/Desktop/biotools/plink/plink --vcf test.vcf --make-bed --out test \
--double-id
```

```python
/Users/ico/Desktop/biotools/plink/plink --bfile test --recode --tab --out test
```

```python
bcftools view -S train_id.txt X.vcf > train.vcf
```

```python
/Users/ico/Desktop/biotools/plink/plink --vcf train.vcf --make-bed \
--out train --double-id
```

```python

commit("test")
commit("train")
```

```bash
/Users/ico/Desktop/biotools/plink/plink --bfile train --recode --tab --out train
```

```bash
/Users/ico/Desktop/biotools/plink/plink --bfile train --allow-no-sex --fisher \
--out train
```

```bash
/Users/ico/Desktop/biotools/gtool -P --ped train.ped --map train.map --og train.gen --os train_fake.sample             
```

```bash
bcftools convert train.vcf -g trains
```

```bash
! /content/snptest_v2.5.2_linux_x86_64_static/snptest_v2.5.2 -data trains.gen.gz \
train_snptest.sample -o train.sum -frequentist 1 -method score -pheno pheno
```

```bash
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
```

```bash
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
```

```bash
gunzip -c Data.txt.gz | awk 'NR==1 || ($11 > 0.01) && ($10 > 0.8) {print}' \
| gzip > Data.gz
```

```bash
gunzip -c Data.gz | awk '{seen[$3]++; if(seen[$3]==1){print}}' \
| gzip -> Data.nodup.gz
```

```bash
gunzip -c Data.nodup.gz |\
awk '!( ($4=="A" && $5=="T") || \
        ($4=="T" && $5=="A") || \
        ($4=="G" && $5=="C") || \
        ($4=="C" && $5=="G")) {print}' |\
    gzip > Data.QC.gz
```

```bash
/Users/ico/Desktop/biotools/plink/plink --bfile test --maf 0.01 --hwe 1e-6 \
--geno 0.01 --mind 0.01 --write-snplist --make-just-fam --allow-no-sex \
-out test.QC
```

```bash
/Users/ico/Desktop/biotools/plink/plink --bfile test --keep test.QC.fam \
--extract test.QC.snplist --allow-no-sex --indep-pairwise 200 50 0.25 \
--out test.QC
```

```bash
/Users/ico/Desktop/biotools/plink/plink --bfile test \
--extract test.QC.prune.in \
--keep test.QC.fam \
--het \
--out test.QC
```

```bash
os.system("Rscript QCtarget.R "+direc+" 1")
print("Rscript QCtarget.R 1")
```

```bash
/Users/ico/Desktop/biotools/plink/plink --bfile test \
--extract test.QC.prune.in \
--make-bed \
--keep test.valid.sample \
--out test.QC
```

```bash
/Users/ico/Desktop/biotools/plink/plink --bfile test \
--make-bed --allow-no-sex \
--out test.QC \
--extract test.QC.snplist \
--exclude test.mismatch \
--a1-allele test.a1
```

```bash
/Users/ico/Desktop/biotools/plink/plink --bfile test \
--make-bed --allow-no-sex --out test.QC \
--extract test.QC.snplist --a1-allele test.a1
```

```bash
os.system("Rscript QCtarget.R "+direc+" 3")
```

```bash
/Users/ico/Desktop/biotools/plink/plink  --bfile test \
--clump-p1 1 --clump-r2 0.1 --clump-kb 250 \
--clump Data.QC.Transformed --clump-snp-field SNP --clump-field P --out test
```

```bash
awk 'NR!=1{print $3}' test.clumped > test.valid.snp
```

```bash
awk '{print $3,$8}' Data.QC.Transformed > SNP.pvalue
```

```bash
echo \"0.001 0 0.001\" > range_list
echo \"0.05 0 0.05\" >> range_list
echo \"0.1 0 0.1\" >> range_list
echo \"0.2 0 0.2\" >> range_list
echo \"0.3 0 0.3\" >> range_list
echo \"0.4 0 0.4\" >> range_list
echo \"0.5 0 0.5\" >> range_list
```

```bash
/Users/ico/Desktop/biotools/plink/plink --bfile test --score Data.QC.Transformed 3 4 9 header --q-score-range range_list SNP.pvalue --extract test.valid.snp --out test
```

```bash
/Users/ico/Desktop/biotools/plink/plink  --bfile test.QC --indep-pairwise 200 50 0.25 --out test
```

```bash
/Users/ico/Desktop/biotools/plink/plink --bfile test.QC --extract test.prune.in --pca 6 --out test
```

```bash
# Read the QC GWAS.
pvalues = pd.read_csv("Data.QC.gz",compression='gzip',sep="\t")
print(pvalues.head())
# The selection of p-value is explained in the manuscript. We started from the lower p-value threshold and then moved to significant p-value thresholds.
pvalue = float(0.05)
pvalues['P']=pd.to_numeric(pvalues['P'],errors='coerce')
subpvalues = pvalues[pvalues['P']<float(0.05)]
subpvalues.to_csv('pv_0.05'+'.txt', columns=['SNP'],index=False,header=False)
```

```bash
/Users/ico/Desktop/biotools/plink/plink --bfile test --extract pv_0.05.txt \
--recodeA --out pv_0.05.ptest
```

```bash
/Users/ico/Desktop/biotools/plink/plink --bfile train --extract pv_0.05.txt \
--recodeA --out pv_0.05.ptrain
```

All previous steps we make for all groups - European, African and East Asian.

Finally, we get necessary files for ML algorithms.

Code for Machine Learning is kinda large, and that’s why only the metrics are presented, but it’s available here:

[Google Colaboratory](https://colab.research.google.com/drive/1LepiYd-p9VvkcG22C3SScCYTrYaQaI2b?usp=sharing)

The metrics **when the model learned on one population predicts another one**:

The first one is the population, on which the model was learned, the second one is the population what the model predicts. 

![Untitled](Simulation%20of%20genotype-phenotype%20for%20different%20pop%20bd860396e0c748fc9389028068de0706/Untitled%205.png)

The metrics when the model is learned on **all populations** and predicts one population:

![Untitled](Simulation%20of%20genotype-phenotype%20for%20different%20pop%20bd860396e0c748fc9389028068de0706/Untitled%206.png)

## SimulatePhenotypes from HapGen

*In process…*

## GEPSI

*In process…*
