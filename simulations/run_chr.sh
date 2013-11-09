#!/bin/bash -vx

#This is the same as run_simple_simulation.sh, except it takes an argument for which chromosome you want 
# to run on and the lengt. so see the comments for that file and run this like ./run_simple_chr.sh 22 60e6 etc... 
# Then, to run on all chromosomes, edit the entries below to point where you want, and run run_all.sh
# which will run this for chrs 1...22 and then run some additional scripts to consolidate the results. 

##############################################################################################################
# Script arguments
CHR=$1
nbp=$2
sim_type=$3
sleep_time=$4

##############################################################################################################

# Edit simulation parameters here. 

# where is macs?
MACS_DIR=~/Packages/macs
# Where do you want the simulations to go?
SIMS_DIR=/data1/users/mathii/f2_sims/${sim_type}/chr${CHR}
# Where is the recombination map, in impute format?
HM2_MAP=~/hm2_recombination_map/genetic_map_GRCh37_chr${CHR}.txt.gz
# Where is the code - this point to the directory you downloaded from github
CODE_DIR=~/f2/f2_age

# Parameters: number of hapotypes, Ne, estimated doubleton power, mutation rate 
nhp1=100
nhp2=100
nhp=200
ne=14000
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

sleep ${sleep_time}

# redirect output to logfile                                                                                                                                   
LOG=${RD}/log.txt
# exec > ${LOG} 2>&1

# compound params
theta=`echo "4*$ne*$mu" | bc`
rho=`echo "4*$ne*0.00000001" | bc`

# 1) convert hms map to macs format
R --vanilla --args ${HM2_MAP} ${MD}/map.txt ${MD}/cut.map.txt < ${CD}/scripts/convert_HM_maps_to_macs_format.R

echo $nhp > $RD/groups.txt
echo "NA" > $RD/events.txt

# 2) Simulate using macs 
case $sim_type in
simple)
	${MACS_DIR}/macs ${nhp} ${nbp} -T -t ${theta} -r ${rho} -h 1e3 -R ${MD}/map.txt 2> \
	    ${SIMS_DIR}/raw_macs_data/trees.txt | ${MACS_DIR}/msformatter | gzip -cf > ${SIMS_DIR}/raw_macs_data/haplotypes.txt.gz
	;;
constant_rate)
	${MACS_DIR}/macs ${nhp} ${nbp} -T -t ${theta} -r ${rho} -h 1e3 -R ${CD}/test/map.txt 2> \
	    ${SIMS_DIR}/raw_macs_data/trees.txt | ${MACS_DIR}/msformatter | gzip -cf > ${SIMS_DIR}/raw_macs_data/haplotypes.txt.gz
	;;
expanding)
        ${MACS_DIR}/macs ${nhp} ${nbp} -t ${theta} -r ${rho} \
            -h 1e3 -R ${MD}/map.txt -eN 0 1 -G 10 -T 2> \
            ${SIMS_DIR}/raw_macs_data/trees.txt | ${MACS_DIR}/msformatter \
            | gzip -cf > ${SIMS_DIR}/raw_macs_data/haplotypes.txt.gz
	;;
ss_test)
	ne=20000
	${MACS_DIR}/macs ${nhp} ${nbp} -t 0.00096 -r 0.0008 -I 2 ${nhp1} ${nhp2} \
	    -ej .0041 2 1 -T 2> \
            ${SIMS_DIR}/raw_macs_data/trees.txt | ${MACS_DIR}/msformatter \
            | gzip -cf > ${SIMS_DIR}/raw_macs_data/haplotypes.txt.gz
;;
bottleneck)
	echo "bottleneck 15" > $RD/events.txt
        ${MACS_DIR}/macs ${nhp} ${nbp} -t ${theta} -r ${rho} \
            -h 1e3 -R ${MD}/map.txt -eN 0.00025 0.01 -eN 0.0003 1 -T 2> \
            ${SIMS_DIR}/raw_macs_data/trees.txt | ${MACS_DIR}/msformatter \
            | gzip -cf > ${SIMS_DIR}/raw_macs_data/haplotypes.txt.gz
	;;
schaffner)
	nhp=400
	ne=100000
	python ${CD}/scripts/schaffner.py $nhp $nbp  \
	    ${MD}/map.txt 12345 ${WD} ${MACS_DIR}
;;
complex)
	nhp=`echo "$nhp1+$nhp2" | bc`
	echo "$nhp1\n$nhp2" > $RD/groups.txt
	echo "split 560" > $RD/events.txt
	${MACS_DIR}/macs ${nhp} ${nbp} -I 2 ${nhp1} ${nhp2} -t ${theta} -r ${rho} \
	    -h 1e3 -R ${MD}/map.txt -n 2 10 -g 2 23.026 -ej 0.01 2 1 -T 2> \
	    ${SIMS_DIR}/raw_macs_data/trees.txt | ${MACS_DIR}/msformatter \
	    | gzip -cf > ${SIMS_DIR}/raw_macs_data/haplotypes.txt.gz
	;;
Ne_10)
	theta10=`echo "4*10*$ne*$mu" | bc`
	rho10=`echo "4*10*$ne*0.00000001" | bc`
        ${MACS_DIR}/macs ${nhp} ${nbp} -T -t ${theta10} -r ${rho10} -h 1e3 -R ${MD}/map.txt 2> \
            ${SIMS_DIR}/raw_macs_data/trees.txt | ${MACS_DIR}/msformatter \
            | gzip -cf > ${SIMS_DIR}/raw_macs_data/haplotypes.txt.gz
        ;;
Ne_01)
        theta01=`echo "4*0.1*$ne*$mu" | bc`
        rho01=`echo "4*0.1*$ne*0.00000001" | bc`
        ${MACS_DIR}/macs ${nhp} ${nbp} -T -t ${theta01} -r ${rho01} -h 1e3 -R ${MD}/map.txt 2> \
            ${SIMS_DIR}/raw_macs_data/trees.txt | ${MACS_DIR}/msformatter \
            | gzip -cf > ${SIMS_DIR}/raw_macs_data/haplotypes.txt.gz
        ;;
ancient_split)
	nhp=`echo "$nhp1+$nhp2" | bc`
	echo "$nhp1\n$nhp2" > $RD/groups.txt
        ${MACS_DIR}/macs ${nhp} ${nbp} -I 2 ${nhp1} ${nhp2} -t ${theta} -r ${rho} \
            -h 1e3 -R ${MD}/map.txt -ej 0.02 2 1 -T 2> \
            ${SIMS_DIR}/raw_macs_data/trees.txt | ${MACS_DIR}/msformatter \
            | gzip -cf > ${SIMS_DIR}/raw_macs_data/haplotypes.txt.gz
	;;
ancient_split_migration)
	nhp=`echo "$nhp1+$nhp2" | bc`
	echo "$nhp1\n$nhp2" > $RD/groups.txt
        ${MACS_DIR}/macs ${nhp} ${nbp} -I 2 ${nhp1} ${nhp2} -t ${theta} -r ${rho} \
            -h 1e3 -R ${MD}/map.txt -ej 0.02 2 1 -ema 0.01 2 x 560 560 x -T 2> \
            ${SIMS_DIR}/raw_macs_data/trees.txt | ${MACS_DIR}/msformatter \
            | gzip -cf > ${SIMS_DIR}/raw_macs_data/haplotypes.txt.gz
	;;
ancient_split_migration_growth)
	nhp=`echo "$nhp1+$nhp2" | bc`
	ne=140000
	theta=`echo "4*$ne*$mu" | bc`
	rho=`echo "4*$ne*0.00000001" | bc`
	echo "$nhp1\n$nhp2" > $RD/groups.txt
        ${MACS_DIR}/macs ${nhp} ${nbp} -I 2 ${nhp1} ${nhp2} -t ${theta} -r ${rho} \
            -h 1e3 -R ${MD}/map.txt  -eN 0.00201 0.1 -ej 0.002 2 1 \
	    -ema 0.001 2 x 560 560 x -G 115 2> \
            ${SIMS_DIR}/raw_macs_data/trees.txt | ${MACS_DIR}/msformatter \
            | gzip -cf > ${SIMS_DIR}/raw_macs_data/haplotypes.txt.gz
	;;
ancient_split_migration_bottleneck)
	nhp=`echo "$nhp1+$nhp2" | bc`
	echo "$nhp1\n$nhp2" > $RD/groups.txt
        ${MACS_DIR}/macs ${nhp} ${nbp} -I 2 ${nhp1} ${nhp2} -t ${theta} -r ${rho} \
            -h 1e3 -R ${MD}/map.txt -ej 0.02 2 1 -ema 0.01 2 x 560 560 x \
	    -en 0.0198 1 1 -en 0.0196 1 0.01 -en 0.0197 2 1 -en 0.0195 2 0.01 -T 2> \
            ${SIMS_DIR}/raw_macs_data/trees.txt | ${MACS_DIR}/msformatter \
            | gzip -cf > ${SIMS_DIR}/raw_macs_data/haplotypes.txt.gz
	;;
ancient_split_long_migration)
	nhp=`echo "$nhp1+$nhp2" | bc`
	echo "$nhp1\n$nhp2" > $RD/groups.txt
        ${MACS_DIR}/macs ${nhp} ${nbp} -I 2 ${nhp1} ${nhp2} -t ${theta} -r ${rho} \
            -h 1e3 -R ${MD}/map.txt -ej 0.02 2 1 -ema 0.004 2 x 560 560 x -T 2> \
            ${SIMS_DIR}/raw_macs_data/trees.txt | ${MACS_DIR}/msformatter \
            | gzip -cf > ${SIMS_DIR}/raw_macs_data/haplotypes.txt.gz
	;;
esac

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
R --vanilla --quiet --slave --args ${CD} ${TH} ${RD} ${TH}/samples.txt ${MD}/cut.map.txt 100 100 two.way ${nbp} < ${CD}/scripts/estimate_error_parameters.R
R --vanilla --quiet --slave --args ${CD} ${RD} ${ne} ${nhp} ${mu} 6 60 < ${CD}/scripts/compare_estimates.R
