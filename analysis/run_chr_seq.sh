#!/bin/bash -vx

# Run analysis of [1000 genomes or other] real data. Find all shared haplotypes
# around f2 variants and estimate their age. For a single chromosome. 

##############################################################################################################

CHR=$1
path_to_seq_data=$2
path_to_chip_data=${path_to_seq_data}
nbp=$3

##############################################################################################################

# Edit parameters here. 

# Where do you want the simulations to go?
VCFTOOLS=~/Packages/vcftools_0.1.12a/bin/vcftools
SIMS_DIR=/groups/reich/iain/f2/1000gP3/chr${CHR}
# Where are the recombination maps, in impute format
HM2_MAP=~/recombination_maps/hm2/genetic_map_GRCh37_chr${CHR}.txt.gz
# Where is the code - this point to the directory you downloaded from github
CODE_DIR=~/f2_age/code
# two colum tab delim file - col1=sample name, col2=population
PANEL=${CODE_DIR}/analysis/1kgP3_panel.txt
# Setup file - see the 1kg file for an example
SETUP_FILE=${CODE_DIR}/analysis/1kgP3_setup.R
# Low complexity regions file
lcr_file=/groups/reich/iain/f2/files/hs37d5-LCRs.20140224.bed

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

# 1) Parse data into correct format 
# 1.1) Extract sample names
cut -f1 ${PANEL} > ${TH}/samples.txt

# 1.2) Chip data
# set max number of file descriptors
ulimit -n 3000
# Keep only the samples we want - can skip this step if you know that you
# are keeping all the samples and sites in the chip data. 
${VCFTOOLS} --gzvcf ${path_to_chip_data} --keep ${TH}/samples.txt \
    --remove-indels --exclude-bed ${lcr_file} \
    --out ${TH}/chr${CHR}.1092.tmp --recode-to-stream | gzip -c \
    > ${TH}/chr${CHR}.1092.tmp.vcf.gz
# Convert to flat format
zgrep -v "^#" ${TH}/chr${CHR}.1092.tmp.vcf.gz \
   | awk '{{printf "%d",$2};for(i=10; i<=NF; i++){printf "\t%d\t%d",substr($i,1,1),substr($i,3,1)};printf "\n"}'  \
   | gzip -c > ${TH}/all.chr${CHR}.chip.haps.gz
# Split up by sample. 
python ${CD}/scripts/haps_to_gt_bysample.py -h ${TH}/all.chr${CHR}.chip.haps.gz \
    -o ${TH}/by_sample/ -s ${TH}/samples.txt
# Cleanup
rm ${TH}/chr${CHR}.1092.tmp.log
rm ${TH}/chr${CHR}.1092.tmp.vcf.gz
rm ${TH}/all.chr${CHR}.chip.haps.gz

# 1.3) Sequence data - extract singletons and doubletons
for n in 1 2
do
   ${VCFTOOLS} --gzvcf ${path_to_seq_data}  --out ${TH}/chr${CHR}.f${n}.log.tmp \
       --remove-indels --recode-to-stream --mac $n --max-mac $n --exclude-bed ${lcr_file} \
       | grep -v "^#" \
       | awk '{{printf "%d",$2};for(i=10; i<=NF; i++){printf "\t%d\t%d",substr($i,1,1),substr($i,3,1)};printf "\n"}'  \
       | gzip -c > ${TH}/chr${CHR}.f${n}.haps.gz

   python ${CD}/scripts/calculate_fn_sites.py -h ${TH}/chr${CHR}.f${n}.haps.gz \
       -o ${TH}/pos.idx.f${n}.gz -n $n > ${TH}/pos.idx.f${n}.tmp.log 
done

# 2) Find f2 haplotypes
R --vanilla --quiet --slave --args ${CD} ${TH} ${RD}/f2_haplotypes.txt \
    ${HM2_MAP} 0 8 < ${CD}/scripts/haplotypes_from_f2.R
gzip -f ${RD}/f2_haplotypes.txt 

# 3) Estimate distributions
R --vanilla --quiet --slave --args ${CD} ${TH} ${RD} ${TH}/samples.txt \
    ${HM2_MAP} 100 100 two.way ${nbp} < ${CD}/scripts/estimate_error_parameters.R
R --vanilla --quiet --slave --args ${CD} ${RD} ${ne} 1092 ${mu} 6 60 \
    ${SETUP_FILE} < ${CD}/scripts/run_analysis_chr.R
