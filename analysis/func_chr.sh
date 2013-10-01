#!/bin/bash -vx

# Run analysis of [1000 genomes or other] real data. Assume we already ran run_chr and have generated
# results for the full dataset. Now we split up by functional class. 

##############################################################################################################

CHR=$1
path_to_seq_data=$2
path_to_chip_data=$3
nbp=$4

##############################################################################################################

# Edit parameters here. 

# Where do you want the simulations to go?
VCFTOOLS=~/Packages/vcftools_0.1.11/bin/vcftools
SIMS_DIR=/data1/users/mathii/1000g/results/chr${CHR}
# Where are the recombination maps, in impute format
HM2_MAP=~/hm2_recombination_map/genetic_map_GRCh37_chr${CHR}.txt.gz
# Where is the code - this point to the directory you downloaded from github
CODE_DIR=~/f2/f2_age

# Parameters: number of hapotypes, Ne, estimated doubleton power, mutation rate 
ne=18000
mu=0.000000004

##############################################################################################################

TH=${SIMS_DIR}/haplotypes
RD=${SIMS_DIR}/results
CD=${CODE_DIR}
PATH_TO_FUNC=/data1/users/mathii/1000g/data/func/ # ALL_XX_pos.txt.gz where XX=coding,noncoding and lof
SAMPLES=${CODE_DIR}/analysis/1092_phase1_samples.txt
for dir in ${SIMS_DIR}/coding ${SIMS_DIR}/noncoding ${SIMS_DIR}/lof
do
    mkdir -p ${dir}/haplotypes
    mkdir -p ${dir}/results
done

LOG=${RD}/log.txt
exec > ${LOG} 2>&1                                                                                                                                        
theta=`echo "4*$ne*$mu" | bc`

for func in lof coding noncoding
do
    echo ${func}
    FTH=${SIMS_DIR}/${func}/haplotypes
    FRD=${SIMS_DIR}/${func}/results
# 1) Sequence data - doubletons from sequences filted by functional class. 
    ${VCFTOOLS} --gzvcf ${path_to_seq_data} --out ${FTH}/chr${CHR}.f2.log.tmp \
        --recode-to-stream --mac 2 --max-mac 2 --positions ${PATH_TO_FUNC}/ALL_${func}_pos.txt.gz | grep -v "^#" \
        | awk '{{printf "%d",$2};for(i=10; i<=NF; i++){printf "\t%d\t%d",substr($i,1,1),substr($i,3,1)};printf "\n"}'  \
        | gzip -c > ${FTH}/chr${CHR}.f2.haps.gz

    python ${CD}/scripts/calculate_fn_sites.py -h ${TH}/chr${CHR}.f2.haps.gz -o ${FTH}/pos.idx.f2.gz -n 2 > ${FTH}/pos.idx.f2.tmp.log 
    ln -sf ${TH}/pos.idx.f1.gz ${FTH}/pos.idx.f1.gz # symlink to singletons
    ln -sf ${TH}/by_sample ${FTH}/by_sample
    ln -sf ${TH}/samples.txt ${FTH}/samples.txt

# 2) Find f2 haplotypes
    R --vanilla --quiet --slave --args ${CD} ${FTH} ${FRD}/f2_haplotypes.txt ${HM2_MAP} 0 < ${CD}/scripts/haplotypes_from_f2.R
    gzip -f ${FRD}/f2_haplotypes.txt 

# 3) Estimate distributions
    ln -sf ${RD}/error_params.txt ${FRD}/error_params.txt
    ln -sf ${RD}/theta_estimates.txt ${FRD}/theta_estimates.txt

    R --vanilla --quiet --slave --args ${CD} ${FRD} ${ne} 1092 ${mu} 6 60 < ${CD}/scripts/run_1kg_analysis_chr.R
done

