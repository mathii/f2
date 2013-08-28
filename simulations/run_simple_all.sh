#!/bin/sh

# Run the script run_simple_chr.sh for all chromosomes, then run some additional scripts to
# consolitdate the results and report on the whole genome results. Need to edit the parameters
# inside run_simple_chr.sh and tell this script where to look for the results. 

##############################################################################################################

# Where are the simulation results. Should be of the form SIMS_ROOT/chr${CHR}/ where chr runs from 1 to 22
SIMS_ROOT=~/f2/simulations/simple/
# Where is the code - this should point to the directory you downloaded from github
CODE_DIR=~/f2/f2_age
# How many parallel processes?
N_PROCS=8
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
    echo "$CODE_DIR/run_simple_chr.sh ${CHR} ${CHR_LENGTHS[${CHR}]}" >> ${SIMS_ROOT}/tmp_args
done

# xargs --arg-file=${SIMS_ROOT}/tmp_args --max-procs=N_PROCS --replace --verbose /bin/sh -c "{}"
# rm ${SIMS_ROOT}/tmp_args

# # Now run files to combine results
# mkdir -p ${SIMS_ROOT}/all

