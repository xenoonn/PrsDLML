cut -f 550-743 -d " " /home/kkoshkina/simulation/hapmap3_r2_b36/hapmap3_r2_b36_chr22.haps > chr22_1.afr_subset.hap
cut -f 898-1012 -d " " /home/kkoshkina/simulation/hapmap3_r2_b36/hapmap3_r2_b36_chr22.haps > chr22_2.afr_subset.hap

cat chr22_1.afr_subset.hap chr22_2.afr_subset.hap > chr22.afr.hap

progdir=/home/kkoshkina/simulation/ 

for chr in $(seq 1 22); do
        dummyDL=`sed -n '2'p /home/kkoshkina/simulation/hapmap3_r2_b36/hapmap3_r2_b36_chr${chr}.legend | cut -d ' ' -f 2`
        $progdir/hapgen2 -m /home/kkoshkina/simulation/hapmap3_r2_b36/genetic_map_chr${chr}_combined_b36.txt \
                -l /home/kkoshkina/simulation/hapmap3_r2_b36/hapmap3_r2_b36_chr${chr}.legend \
                -h /home/kkoshkina/simulation/eur/chr${chr}.afr.hap \
                -o /home/kkoshkina/simulation/eur/genotypes_chr${chr}_hapgen \
                -n 5000 0 \
                -dl $dummyDL 0 0 0 \
                -no_haps_output
done
