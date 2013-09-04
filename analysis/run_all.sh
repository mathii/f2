#!/bin/bash

# Run the 1000 genomes analysis for all chromosomes and then combine the results

##############################################################################################################

# Where are the simulation results. Should be of the form SIMS_ROOT/chr${CHR}/ where chr runs from 1 to 22
SIMS_ROOT=~/f2/simulations/analysis/
# Where is the code - this should point to the directory you downloaded from github
CODE_DIR=~/f2/f2_age
# How many parallel processes?
N_PROCS=22

##############################################################################################################

# from http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/data/
# First element is dummy... 
declare -a CHR_LENGTHS=(-1 249250621 243199373 198022430 191154276 180915260 171115067\
 159138663 146364022 141213431 135534747 135006516 133851895 115169878 107349540\
 102531392 90354753 81195210 78077248 59128983 63025520 48129895 51304566) 

##############################################################################################################

mkdir -p ${SIMS_ROOT}

rm -f ${SIMS_ROOT}/tmp_args
for CHR in {1..22}
do
    echo "$CODE_DIR/analysis/run_chr.sh ${CHR} /data1/users/mathii/1000g/data/seq/ALL.chr${CHR}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz /data1/users/mathii/1000g/data/chip/ALL.chr${CHR}.omni_2123_samples_b37_SHAPEIT.20120103.snps.chip_based.haplotypes.vcf.gz ${CHR_LENGTHS[${CHR}]}" >> ${SIMS_ROOT}/tmp_args
done

xargs --arg-file=${SIMS_ROOT}/tmp_args --max-procs=${N_PROCS} --replace --verbose /bin/bash -c "{}"
rm ${SIMS_ROOT}/tmp_args

# Now run files to combine results
mkdir -p ${SIMS_ROOT}/all


 