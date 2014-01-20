#!/bin/bash -vx

# Run analysis of [1000 genomes or other] real data. Find all shared haplotypes
# around fn for n=3 to 5 variants and estimate their age. For a single chromosome. 
# Assume 2 has already been done. 				

##############################################################################################################

CHR=$1
path_to_seq_data=$2

##############################################################################################################

# Edit parameters here. 

# Where do you want the simulations to go?
VCFTOOLS=~/Packages/vcftools_0.1.11/bin/vcftools
SIMS_DIR=/data1/users/mathii/1000g/results/chr${CHR}
# Where are the recombination maps, in impute format
HM2_MAP=~/hm2_recombination_map/genetic_map_GRCh37_chr${CHR}.txt.gz
# Where is the code - this point to the directory you downloaded from github
CODE_DIR=~/f2/f2_age
# two colum tab delim file - col1=sample name, col2=population
PANEL=${CODE_DIR}/analysis/1kg_panel.txt
# Setup file - see the 1kg file for an example
SETUP_FILE=${CODE_DIR}/analysis/1kg_setup.R

# Parameters: number of hapotypes, Ne, estimated doubleton power, mutation rate 
ne=185000
mu=0.000000004

##############################################################################################################

TH=${SIMS_DIR}/haplotypes
RD=${SIMS_DIR}/results
CD=${CODE_DIR}

for dir in ${WD} ${TH}/by_sample ${RD}
do
    mkdir -p ${dir}
done

LOG=${RD}/log.txt
exec > ${LOG} 2>&1                                                                                                                                        
# assume that run_chr has already been run. 

# 1.3) Sequence data - extract variants, find and date haplotypes
for n in 3 4 5
do
   ${VCFTOOLS} --gzvcf ${path_to_seq_data}  --out ${TH}/chr${CHR}.f${n}.log.tmp \
       --recode-to-stream --mac $n --max-mac $n | grep -v "^#" \
       | awk '{{printf "%d",$2};for(i=10; i<=NF; i++){printf "\t%d\t%d",substr($i,1,1),substr($i,3,1)};printf "\n"}'  \
       | gzip -c > ${TH}/chr${CHR}.f${n}.haps.gz

   python ${CD}/scripts/calculate_fn_sites.py -h ${TH}/chr${CHR}.f${n}.haps.gz \
       -o ${TH}/pos.idx.f${n}.gz -n $n > ${TH}/pos.idx.f${n}.tmp.log 

   R --vanilla --quiet --slave --args ${CD} ${TH} ${RD} \
       ${HM2_MAP} 0 ${n} ${ne} < ${CD}/scripts/haplotypes_from_fn.R
   gzip -f ${RD}/f${n}_results.txt 
done

