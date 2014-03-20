# Generate macs data, of with specified parameters. 
run_sim_type=$3

case $run_sim_type in
    simple)
	${MACS_DIR}/macs ${nhp} ${nbp} -T -t ${theta} -r ${rho} -h 1e3 -R ${MD}/map.txt 2> \
	    ${SIMS_DIR}/raw_macs_data/trees.txt | ${MACS_DIR}/msformatter | gzip -cf > ${SIMS_DIR}/raw_macs_data/haplotypes.txt.gz
	;;
constant_rate)
	${MACS_DIR}/macs ${nhp} ${nbp} -T -t ${theta} -r ${rho} -h 1e3 -R ${MD}/map.txt 2> \
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
bottleneck_10)
        ${MACS_DIR}/macs ${nhp} ${nbp} -t ${theta} -r ${rho} \
            -h 1e3 -R ${MD}/map.txt -eN 0.0001785714 0.05 -eN 0.0002678571 1 -T 2> \
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
split_10)
	# time = 
	nhp=`echo "$nhp1+$nhp2" | bc`
	echo "$nhp1\n$nhp2" > $RD/groups.txt
        ${MACS_DIR}/macs ${nhp} ${nbp} -I 2 ${nhp1} ${nhp2} -t ${theta} -r ${rho} \
            -h 1e3 -R ${MD}/map.txt -ej 0.0001785714 2 1 -T 2> \
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
ancient_split_wrong_map)
	R --vanilla --args "~/aa_recombination_map/maps_chr_hg19.${CHR}.gz" ${MD}/map.txt \
	    < ${CD}/scripts/convert_AA_maps_to_macs_format.R
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
ancient_split_growth)
	nhp=`echo "$nhp1+$nhp2" | bc`
	ne=140000
	theta=`echo "4*$ne*$mu" | bc`
	rho=`echo "4*$ne*0.00000001" | bc`
	echo "$nhp1\n$nhp2" > $RD/groups.txt
        ${MACS_DIR}/macs ${nhp} ${nbp} -I 2 ${nhp1} ${nhp2} -t ${theta} -r ${rho} \
            -h 1e3 -R ${MD}/map.txt  -eN 0.00201 0.1 -ej 0.002 2 1 -G 115 -T 2> \
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
	    -ema 0.001 2 x 560 560 x -G 115 -T 2> \
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
ancient_split_migration_interval)
	nhp=`echo "$nhp1+$nhp2" | bc`
	echo "$nhp1\n$nhp2" > $RD/groups.txt
        ${MACS_DIR}/macs ${nhp} ${nbp} -I 2 ${nhp1} ${nhp2} -t ${theta} -r ${rho} \
            -h 1e3 -R ${MD}/map.txt -ej 0.02 2 1 -ema 0.01 2 x 0 0 x \
	     -ema 0.005 2 x 560 560 x 2> \
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
