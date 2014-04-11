#!/bin/bash -vx

# Run analysis of [1000 genomes or other] real data. Find all shared haplotypes
# around fn for n=3 to 5 variants and estimate their age. For a single chromosome. 
# Assume 2 has already been done. 				

##############################################################################################################

CHR=$1
nbp=$2
sim_type=$3

##############################################################################################################

# Edit parameters here. 

# Where do you want the simulations to go?
SIMS_DIR=/data/mathii/f2/simulations/${sim_type}/chr${CHR}
# # Where are the recombination maps, in impute format
HM2_MAP=~/recombination_maps/hm2/genetic_map_GRCh37_chr${CHR}.txt.gz
# Where is the code - this point to the directory you downloaded from github
CODE_DIR=~/f2/f2_age

# Parameters: number of hapotypes, Ne, estimated doubleton power, mutation rate 
nhp1=200
nhp2=200
nhp=200
ne=14000
mu=0.000000012

##############################################################################################################

WD=${SIMS_DIR}/raw_macs_data
TH=${SIMS_DIR}/haplotypes
RD=${SIMS_DIR}/results
CD=${CODE_DIR}
                           
# assume that run_chr has already been run. 

# 1.3) Sequence data - extract variants, find and date haplotypes
# for comparison, do n=2 in the same way - sanity, check it agrees
# with the full f2 analysis
n=2
R --vanilla --quiet --slave --args ${CD} ${TH} ${RD} \
    ${HM2_MAP} 1 ${n} ${ne} < ${CD}/scripts/haplotypes_from_fn.R
gzip -f ${RD}/f${n}_results.txt 

# Now extract haplotypes for f3 to f5
python ${CD}/scripts/macs_genotype_to_hap_files.py -g ${WD}/genotypes.txt.gz \
    -p ${WD}/snps.pos.txt.gz -l ${nbp} -o ${TH} -f 3 -t 5

for n in 3 4 5
do
   python ${CD}/scripts/calculate_fn_sites.py -h ${TH}/haps.f${n}.gz \
       -o ${TH}/pos.idx.f${n}.gz -n $n > ${TH}/pos.idx.f${n}.tmp.log 

   R --vanilla --quiet --slave --args ${CD} ${TH} ${RD} \
       ${HM2_MAP} 1 ${n} ${ne} < ${CD}/scripts/haplotypes_from_fn.R
   gzip -f ${RD}/f${n}_results.txt 
done

