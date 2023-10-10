
# convert to plink format and prune SNPs
for chr in `seq 1 22`; do
	plink --data genotypes_chr${chr}_hapgen.controls \
		--oxford-single-chr $chr \
		--make-bed \
		--out genotypes_chr${chr}_hapgen.controls


	plink --bfile genotypes_chr${chr}_hapgen.controls \
		--indep-pairwise 50kb 1 0.5 \
		--out genotypes_chr${chr}_hapgen.controls

	plink --bfile genotypes_chr${chr}_hapgen.controls \
		--extract genotypes_chr${chr}_hapgen.controls.prune.in \
		--make-bed \
		--out genotypes_chr${chr}_hapgen.controls.pruned
	echo "genotypes_chr${chr}_hapgen.controls.pruned" >> file_list
done

# Merge chromsome-wide files into a single, genome-wide file
plink --merge-list file_list --make-bed --allow-no-sex --out genotypes_genome_hapgen.controls

# compute kinship
plink --bfile genotypes_genome_hapgen.controls\
	--make-rel square \
	--out genotypes_genome_hapgen.controls.grm

# need to rename to oxford format for Phenotype Simulator
for i in $(seq 1 22); do
	cp genotypes_chr${i}_hapgen.controls.sample genotypes_chr${i}_hapsample.controls.sample 
done