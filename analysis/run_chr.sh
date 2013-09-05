#!/bin/bash

# Run analysis of [1000 genomes or other] real data. Find all shared haplotypes
# around f2 variants and estimate their age. For a single chromosome. 

##############################################################################################################

CHR=$1
path_to_seq_data=$2
path_to_chip_data=$3
nbp=$4

##############################################################################################################

# Edit parameters here. 

# Where do you want the simulations to go?
SIMS_DIR=~/f2/1000g/results/chr${CHR}
# Where are the recombination maps, in impute format
HM2_MAP=~/hm2_recombination_map/genetic_map_GRCh37_chr${CHR}.txt.gz
# Where is the code - this point to the directory you downloaded from github
CODE_DIR=~/f2/f2_age

# Parameters: number of hapotypes, Ne, estimated doubleton power, mutation rate 
ne=18000
mu=0.000000012

##############################################################################################################

TH=${SIMS_DIR}/haplotypes
RD=${SIMS_DIR}/results
CD=${CODE_DIR}
SAMPLES=${CODE_DIR}/analysis/1092_phase1_samples.txt
for dir in ${WD} ${TH}/by_sample ${RD}
do
    mkdir -p ${dir}
done

theta=`echo "4*$ne*$mu" | bc`

# 1) Parse data into correct format 
# 1.1) Chip data
# set max number of file descriptors
ulimit -n 1200
# Keep only the Phase 1 samples
vcftools --gzvcf ${path_to_chip_data} --keep ${SAMPLES} \
    --out ${TH}/chr${CHR}.1092.tmp --recode-to-stream | gzip -c > ${TH}/chr${CHR}.1092.tmp.vcf.gz
# Convert to flat format
zgrep -v "^#" ${TH}/chr${CHR}.1092.tmp.vcf.gz \
    | awk '{{printf "%d",$2};for(i=10; i<=NF; i++){printf "\t%d\t%d",substr($i,1,1),substr($i,3,1)};printf "\n"}'  \
    | gzip -c > ${TH}/all.chr${CHR}.chip.haps.gz
# Split up by sample. 
python ${CD}/scripts/haps_to_gt_bysample.py -h ${TH}/all.chr${CHR}.chip.haps.gz -o ${TH}/by_sample/ -s ${SAMPLES}
# Cleanup
rm ${TH}/chr${CHR}.1092.tmp.log
rm ${TH}/chr${CHR}.1092.tmp.vcf.gz
rm ${TH}/all.chr${CHR}.chip.haps.gz

# 1.2) Sequence data - extract singletons and doubletons
for n in 1 2
do
    vcftools --gzvcf ${path_to_seq_data}  --out ${TH}/chr${CHR}.f${n}.log.tmp \
        --recode-to-stream --mac $n --max-mac $n | grep -v "^#" \
        | awk '{{printf "%d",$2};for(i=10; i<=NF; i++){printf "\t%d\t%d",substr($i,1,1),substr($i,3,1)};printf "\n"}'  \
        | gzip -c > ${TH}/chr${CHR}.f${n}.haps.gz

    python ~/trees/code/scripts/calculate_fn_sites.py -h ${TH}/chr${CHR}.f${n}.haps.gz -o ${TH}/pos.idx.f${n}.gz -n $n > ${TH}/pos.idx.f${n}.tmp.log 
done


# 2) Find f2 haplotypes
cp ${CD}/analysis/1092_phase1_samples.txt ${TH}/samples.txt
R --vanilla --quiet --slave --args ${CD} ${TH} ${RD}/f2_haplotypes.txt ${HM2_MAP} 0 < ${CD}/scripts/haplotypes_from_f2.R
gzip -f ${RD}/f2_haplotypes.txt 

# 3) Estimate distributions
R --vanilla --quiet --slave --args ${CD} ${TH} ${RD} ${SAMPLES} ${HM2_MAP} 20 100 two.way ${nbp} < ${CD}/scripts/estimate_error_parameters.R
R --vanilla --quiet --slave --args ${CD} ${RD} ${ne} 1092 ${mu} 6 60 < ${CD}/scripts/run_1kg_analysis_chr.R