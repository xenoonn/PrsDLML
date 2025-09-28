library(PhenotypeSimulator)
indir <- "{population_dir}"

datadir <- "{population_dir}/phensim"
if (!dir.exists(datadir)) dir.create(datadir, recursive=TRUE)
# specify filenames and parameters
totalGeneticVar <- 0.6
totalSNPeffect <- 0.1
h2s <- totalSNPeffect/totalGeneticVar
kinshipfile <- paste(indir, "/genotypes_genome_hapgen.controls.grm.rel",
                     sep="")
genoFilePrefix <- paste(indir, "/genotypes_", sep="")
genoFileSuffix <- "_hapgen.controls.gen"
# simulate phenotype with one phenotype components
simulation <- runSimulation(N = 5000, P = 1, cNrSNP=30000, seed=43,
                            format = "oxgen", # = TRUE,
                            genoFilePrefix = genoFilePrefix,
                            genoFileSuffix = genoFileSuffix,
                            chr = 1:22,
                            mBetaGenetic = 0, sdBetaGenetic = 0.2,
                            theta=1,
                            genVar = totalGeneticVar, h2s = h2s,
                            phi = 0.6, delta = 0.2, rho=0.2,
                            NrFixedEffects = 2, NrConfounders = c(2, 2),
                            distConfounders = c("bin", "norm"),
                            probConfounders = 0.2,
                            genoDelimiter=" ",
                            kinshipDelimiter="\t",
                            kinshipHeader=FALSE,
                            verbose = TRUE )

outdirectory <- savePheno(simulation, directory = datadir,
                          format=c("csv", "plink"),
                          saveIntermediate=TRUE)
