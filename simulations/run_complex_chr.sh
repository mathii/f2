#!/bin/bash -vx

#This is the same as run_simple_chr.sh, except it runs a more complex simulation with a population split
# followed by expansion, and then runs an additional script to analyse the difference between the two populations. 


##############################################################################################################
# script arguments
CHR=$1
nbp=$2

##############################################################################################################

# Edit simulation parameters here. 

# where is macs?
MACS_DIR=~/Packages/macs
# Where do you want the simulations to go?
SIMS_DIR=~/f2/simulations/complex/chr${CHR}
# Where is the recombination map, in impute format?
HM2_MAP=~/hm2_recombination_map/genetic_map_GRCh37_chr${CHR}.txt.gz
# Where is the code - this point to the directory you downloaded from github
CODE_DIR=~/f2/f2_age

# Parameters: number of hapotypes, Ne, estimated doubleton power, mutation rate 
nhp1=200
nhp2=200
ne=10000
dbp=0.66
mu=0.000000012

##############################################################################################################

# Setup directories
WD=${SIMS_DIR}/raw_macs_data
TH=${SIMS_DIR}/haplotypes
RD=${SIMS_DIR}/results
MD=${SIMS_DIR}/map
CD=${CODE_DIR}

for dir in ${WD} ${TH}/by_sample ${RD} ${MD}
do
    mkdir -p ${dir}
done

# redirect output to logfile 
LOG=${RD}/log.txt
# exec > ${LOG} 2>&1

# compound params
theta=`echo "4*$ne*$mu" | bc`
rho=`echo "4*$ne*0.00000001" | bc`
nhp=`echo "$nhp1+$nhp2" | bc`

# 1) convert hms map to macs format
R --vanilla --args ${HM2_MAP} ${MD}/map.txt ${MD}/cut.map.txt < ${CD}/scripts/convert_HM_maps_to_macs_format.R

# 2) Simulate using macs 
${MACS_DIR}/macs ${nhp} ${nbp} -I 2 ${nhp1} ${nhp2} -t ${theta} -r ${rho} -h 1e3 -R ${MD}/map.txt -n 2 10 -g 2 23.026 -ej 0.1 2 1 -T 2> \
    ${SIMS_DIR}/raw_macs_data/trees.txt | ${MACS_DIR}/msformatter | gzip -cf > ${SIMS_DIR}/raw_macs_data/haplotypes.txt.gz

# 3)Parse the macs output to get trees and genotypes
gzip -f ${WD}/trees.txt
# extract trees and cut out the inital first brackets so that we can use the biopython parser. 
zgrep "^\[" ${WD}/haplotypes.txt.gz | sed -e 's/\[[^][]*\]//g' | gzip -cf > ${WD}/trees.newick.txt.gz
# also store tract lengths - use these to recover positions. 
zgrep -o '\[[0-9]*\]' ${WD}/haplotypes.txt.gz | sed 's/\[//g;s/\]//g' | gzip -cf > ${WD}/trees.lengths.txt.gz
# # Get tree pos positions
# Get genotypes ans positions. 
gunzip -c ${WD}/haplotypes.txt.gz | tail -n ${nhp} | gzip -c > ${WD}/genotypes.txt.gz
zgrep "positions" ${WD}/haplotypes.txt.gz | cut -d " " -f2- | gzip -cf > ${WD}/snps.pos.txt.gz 
# Parse the genotype data into the right format, by haplotypes etc... 
python ${CD}/scripts/macs_genotype_to_hap_files.py -g ${WD}/genotypes.txt.gz \
 -p ${WD}/snps.pos.txt.gz -l ${nbp} -o ${TH}

# 4) Get the real haplotyes
python ${CD}/scripts/find_real_nns.py -t ${WD}/trees.newick.txt.gz  \
    -l ${WD}/trees.lengths.txt.gz -o ${TH}/NN_haplotypes.txt.gz -r 1000

# 5) extract f1 and f2 variants. 
python ${CD}/scripts/calculate_fn_sites.py -h ${TH}/haps.f1.gz -o ${TH}/pos.idx.f1.gz -n 1 > ${TH}/f1.pos.tmp.log
python ${CD}/scripts/calculate_fn_sites.py -h ${TH}/haps.f2.gz -o ${TH}/pos.idx.f2.gz -n 2 > ${TH}/f2.pos.tmp.log
python ${CD}/scripts/haps_to_gt_bysample.py -h ${TH}/haps.gz  -o ${TH}/by_sample/ -s ${TH}/samples.txt

# 6) Estimate haplotypes from f2 variants
R --vanilla --quiet --slave --args ${CD} ${TH} ${RD}/f2_haplotypes.txt ${MD}/cut.map.txt 0 < ${CD}/scripts/haplotypes_from_f2.R
gzip -f ${RD}/f2_haplotypes.txt 

# 7) Compare haplotpyes and compute power, then compare estimates of time. 
R --vanilla --quiet --slave --args ${CD} ${TH} ${RD} ${ne} ${dbp} ${MD}/cut.map.txt ${nhp} < ${CD}/scripts/compare_haplotypes.R
gzip -f ${RD}/matched_haps.txt
R --vanilla --quiet --slave --args ${CD} ${TH} ${RD} ${TH}/samples.txt ${MD}/cut.map.txt 20 100 < ${CD}/scripts/estimate_error_parameters.R
R --vanilla --quiet --slave --args ${CD} ${RD} ${ne} ${nhp} ${mu} 6 60 < ${CD}/scripts/compare_estimates.R
