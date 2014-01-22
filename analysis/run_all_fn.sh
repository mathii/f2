# Probably don't need to run these in parallel just yet, but why not?

SIMS_ROOT=/data1/users/mathii/1000g/results/
SEQ_DATA=/data1/users/mathii/1000g/data/seq/ALL.chrXX.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
N_PROCS=22

for chr in {1..22}
do
    seq_path=${SEQ_DATA/XX/$chr}
    echo "$CODE_DIR/analysis/run_chr_fn.sh ${chr} ${seq_path}" >> ${SIMS_ROOT}/tmp_args_fn
done

xargs --arg-file=${SIMS_ROOT}/tmp_args_fn --max-procs=${N_PROCS} --replace --verbose /bin/bash -c "{}"

