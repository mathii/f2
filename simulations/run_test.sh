#!/bin/bash -vx

#Run the whole simulation and analysis pipeline. This runs and analyses a simulation for a population
#with constant size, using the parameters specified below. Suppplied parameters use the recombination 
# map from chr1 (assuming you downloaded it), and simulates 250mb for 100 individuals (200 chromosomes)
# with constant Ne of 10^4. If you change Ne, you'll need to change the population scaled mutation and
# recombination rates in the macs command line.

##############################################################################################################

# Edit simulation parameters here. 

# where is macs?
MACS_DIR=~/Packages/macs
# Where do you want the simulations to go?
SIMS_DIR=~/f2/simulations/test
# Where is the code - this point to the directory you downloaded from github
CODE_DIR=~/f2/f2_age

# Parameters: size of region,number of hapotypes, Ne, estimated power 
nbp=1000000
nhp=100
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
# Where is the recombination map, in impute format?
HM2_MAP=${CD}/test/constant_map.txt.gz

for dir in ${WD} ${TH}/by_sample ${RD} ${MD}
do
    mkdir -p ${dir}
done


# compound params
theta=`echo "4*$ne*$mu" | bc`
rho=`echo "4*$ne*0.00000001" | bc`

# 1) convert hms map to macs format
R --vanilla --args ${HM2_MAP} ${MD}/map.txt ${MD}/cut.map.txt < ${CD}/scripts/convert_HM_maps_to_macs_format.R

# 2) Simulate using macs 
${MACS_DIR}/macs ${nhp} ${nbp} -T -t ${theta} -r ${rho} -h 1e3 -R ${MD}/map.txt 2> \
    ${SIMS_DIR}/raw_macs_data/trees.txt | ${MACS_DIR}/msformatter | gzip -c > ${SIMS_DIR}/raw_macs_data/haplotypes.txt.gz

# 3)Parse the macs output to get trees and genotypes
gzip -f ${WD}/trees.txt
# extract trees and cut out the inital first brackets so that we can use the biopython parser. 
zgrep "^\[" ${WD}/haplotypes.txt.gz | sed -e 's/\[[^][]*\]//g' | gzip -cf > ${WD}/trees.newick.txt.gz
# also store tract lengths - use these to recover positions. 
zgrep -o '\[[0-9]*\]' ${WD}/haplotypes.txt.gz | sed 's/\[//g;s/\]//g' | gzip -cf > ${WD}/trees.lengths.txt.gz
# # Get tree pos positions
# Get genotypes ans positions. 
gunzip -c ${WD}/haplotypes.txt.gz | tail -n ${nhp} | gzip -cf > ${WD}/genotypes.txt.gz
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
R --vanilla --quiet --slave --args ${CD} ${TH} ${RD}/f2_haplotypes.txt ${MD}/cut.map.txt 1 < ${CD}/scripts/haplotypes_from_f2.R
gzip -f ${RD}/f2_haplotypes.txt 

# 7) Compare haplotpyes and compute power, then compare estimates of time. 
R --vanilla --quiet --slave --args ${CD} ${TH} ${RD} ${ne} ${dbp} ${MD}/cut.map.txt ${nhp} < ${CD}/scripts/compare_haplotypes.R
gzip ${RD}/matched_haps.txt
R --vanilla --quiet --slave --args ${CD} ${TH} ${RD} ${TH}/samples.txt ${MD}/cut.map.txt 20 100 two.way < ${CD}/scripts/estimate_error_parameters.R
R --vanilla --quiet --slave --args ${CD} ${RD} ${ne} ${nhp} ${mu} 6 60 < ${CD}/scripts/compare_estimates.R
