#!/bin/bash
#
#SBATCH --job-name=incr
#SBATCH --output=slurm-incr-%j.o.out
#SBATCH --error=slurm-incr-%j.e.out
#
#SBATCH --ntasks=1
#
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4

###########################################################################################################################################################################################################################
# This script is for illustration purposes only. It should run well with few modifications, mostly having to do with paths
# DONT EXPECT THIS TO RUN WITHOUT ANY CHANGES
###########################################################################################################################################################################################################################


## activate a virtual env
source /path/to/virtualenv/bin/activate
## add some paths to PATH
source ~/.fix_paths

macs2="/absolute/path/to/macs2"
py3='/absolute/path/to/python3'
idr='/absolute/path/to/idr'


from=${SLURM_SUBMIT_DIR}
org=$1
stage=$2
gsize=$3
mask=$4


# The nucleosome free reads have been precomputed and we can just grab them from there
ocd1=/path/to/our_nicely_tidy_data/${org}_${stage}_1_clean.bed
ocd2=/path/to/our_nicely_tidy_data/${org}_${stage}_2_clean.bed

# this is part of the queue job stuff
TMPDIR='/scratch/'${SLURM_JOB_ID}'/'
mkdir $TMPDIR
cd $TMPDIR
mkdir ${TMPDIR}/results
mkdir ${TMPDIR}/data

shuf $ocd1 > ./data/inp1.bed
shuf $ocd2 > ./data/inp2.bed


tr1_lines=$( cat ./data/inp1.bed | wc -l )
tr2_lines=$( cat ./data/inp2.bed | wc -l )

# get the minimum read count between the two replicates
minlines=$(echo "${tr1_lines} ${tr2_lines}" | python -c "import sys; print  min([int(x) for x in next(sys.stdin).split()]) ")


# We effectively redo the entire peak calling pipeline using a subsample every time:
turn=0
for (( samplesize=200000; samplesize<=${minlines}; samplesize+=200000 ))
do 
    current_wd=${TMPDIR}/results/loop_${turn}
    mkdir ${current_wd}

	head -n ${samplesize} ./data/inp1.bed | shuf - > ./data/tr1_clean.bed
	head -n ${samplesize} ./data/inp2.bed | shuf - > ./data/tr2_clean.bed

	# count reads in the two replicates 
	tr1_lines=${samplesize}
	tr2_lines=${samplesize}
	# we'll get half for each pseudoreplicate
	tr1_sample_size=$(( (tr1_lines + 1) / 2 ))
	tr2_sample_size=$(( (tr2_lines + 1) / 2 ))
	#Make the replicate 1 pseudoreplicates
	head -n ${tr1_sample_size} ./data/tr1_clean.bed > ./data/tr1_psr1.bed
	tail -n ${tr1_sample_size} ./data/tr1_clean.bed > ./data/tr1_psr2.bed
	#Make the replicate 2 pseudoreplicates
	head -n ${tr2_sample_size} ./data/tr2_clean.bed > ./data/tr2_psr1.bed
	tail -n ${tr2_sample_size} ./data/tr2_clean.bed > ./data/tr2_psr2.bed
	echo 'Got pseudoreplicates'
	# The Pooled pseudoreplicates
	pooled_lines=$(( tr1_lines + tr2_lines ))
	pool_sample_size=$(( (pooled_lines + 1) / 2 ))
	cat ./data/tr1_clean.bed ./data/tr2_clean.bed | shuf > ./data/pooled.bed
	head -n ${pool_sample_size} ./data/pooled.bed > ./data/pool_psr1.bed
	tail -n ${pool_sample_size} ./data/pooled.bed > ./data/pool_psr2.bed
	echo 'oke, got the pooled pseudoreplicates, ready to call peaks'


	# Peak calling for replicate1, replicate 2, and pooled reads
	${macs2} callpeak --nomodel --keep-dup 'auto' --llocal 10000 --extsize 100 --shift -50  -p 0.015 -f BED -g ${gsize} -t ./data/tr1_clean.bed -n ${org}_${stage}_rep1_nf --outdir ${current_wd} &
	${macs2} callpeak --nomodel --keep-dup 'auto' --llocal 10000 --extsize 100 --shift -50  -p 0.015 -f BED -g ${gsize} -t ./data/tr2_clean.bed -n ${org}_${stage}_rep2_nf --outdir ${current_wd} &
	${macs2} callpeak --nomodel --keep-dup 'auto' --llocal 10000 --extsize 100 --shift -50  -p 0.015 -f BED -g ${gsize} -t ./data/tr1_clean.bed ./data/tr2_clean.bed -n ${org}_${stage}_comb_nf --outdir ${current_wd} &
	wait %1 %2 %3
	echo 'Done with proper replicate peak calling'

	#Peak call for the replicate 1 pseudoreplicates
	${macs2} callpeak --nomodel --keep-dup 'auto' --llocal 10000 --extsize 100 --shift -50  -p 0.015 -f BED -g ${gsize} -t ./data/tr1_psr1.bed -n ${org}_${stage}_rep1_psr1 --outdir ${current_wd} &
	${macs2} callpeak --nomodel --keep-dup 'auto' --llocal 10000 --extsize 100 --shift -50  -p 0.015 -f BED -g ${gsize} -t ./data/tr1_psr2.bed -n ${org}_${stage}_rep1_psr2 --outdir ${current_wd} &
	#Peak call for the replicate 2 pseudoreplicates
	${macs2} callpeak --nomodel --keep-dup 'auto' --llocal 10000 --extsize 100 --shift -50  -p 0.015 -f BED -g ${gsize} -t ./data/tr2_psr1.bed -n ${org}_${stage}_rep2_psr1 --outdir ${current_wd} &
	${macs2} callpeak --nomodel --keep-dup 'auto' --llocal 10000 --extsize 100 --shift -50  -p 0.015 -f BED -g ${gsize} -t ./data/tr2_psr2.bed -n ${org}_${stage}_rep2_psr2 --outdir ${current_wd} &
	wait %1 %2 %3 %4
	echo 'Done with pseudo-replicate peak calling'

	#Peak calling for the pooled pseudoreplicates
	${macs2} callpeak --nomodel --keep-dup 'auto' --llocal 10000 --extsize 100 --shift -50  -p 0.015 -f BED -g ${gsize} -t ./data/pool_psr1.bed -n ${org}_${stage}_pool_psr1 --outdir ${current_wd} &
	${macs2} callpeak --nomodel --keep-dup 'auto' --llocal 10000 --extsize 100 --shift -50  -p 0.015 -f BED -g ${gsize} -t ./data/pool_psr2.bed -n ${org}_${stage}_pool_psr2 --outdir ${current_wd} &
	wait %1 %2
	echo 'Done with pooled replicate peak calling'

	# define some names for ease
	peaks1=${current_wd}/${org}_${stage}_rep1_nf_peaks.narrowPeak
	peaks2=${current_wd}/${org}_${stage}_rep2_nf_peaks.narrowPeak
	peaks3=${current_wd}/${org}_${stage}_comb_nf_peaks.narrowPeak
	peaks_ps1=${current_wd}/${org}_${stage}_rep1_psr1_peaks.narrowPeak
	peaks_ps2=${current_wd}/${org}_${stage}_rep1_psr2_peaks.narrowPeak
	peaks_ps3=${current_wd}/${org}_${stage}_rep2_psr1_peaks.narrowPeak
	peaks_ps4=${current_wd}/${org}_${stage}_rep2_psr2_peaks.narrowPeak
	peaks_pool_ps1=${current_wd}/${org}_${stage}_pool_psr1_peaks.narrowPeak
	peaks_pool_ps2=${current_wd}/${org}_${stage}_pool_psr2_peaks.narrowPeak

	# Clean with rep mask
	bedtools intersect -v -a ${peaks1} -b ${mask} -f 0.33 > ${peaks1}_masked &
	bedtools intersect -v -a ${peaks2} -b ${mask} -f 0.33 > ${peaks2}_masked &
	bedtools intersect -v -a ${peaks3} -b ${mask} -f 0.33 > ${peaks3}_masked &
	wait %1 %2 %3
	bedtools intersect -v -a ${peaks_ps1} -b ${mask} -f 0.33 > ${peaks_ps1}_masked &
	bedtools intersect -v -a ${peaks_ps2} -b ${mask} -f 0.33 > ${peaks_ps2}_masked &
	bedtools intersect -v -a ${peaks_ps3} -b ${mask} -f 0.33 > ${peaks_ps3}_masked &
	bedtools intersect -v -a ${peaks_ps4} -b ${mask} -f 0.33 > ${peaks_ps4}_masked &
	wait %1 %2 %3 %4
	bedtools intersect -v -a ${peaks_pool_ps1} -b ${mask} -f 0.33 > ${peaks_pool_ps1}_masked &
	bedtools intersect -v -a ${peaks_pool_ps2} -b ${mask} -f 0.33 > ${peaks_pool_ps2}_masked &
	wait %1 %2


	sort -k8,8nr ${peaks1}_masked | head -n 500000 > ${peaks1}_masked_f &
	sort -k8,8nr ${peaks2}_masked | head -n 500000 > ${peaks2}_masked_f &
	sort -k8,8nr ${peaks3}_masked | head -n 500000 > ${peaks3}_masked_f &
	wait %1 %2 %3
	sort -k8,8nr ${peaks_ps1}_masked | head -n 500000 > ${peaks_ps1}_masked_f &
	sort -k8,8nr ${peaks_ps2}_masked | head -n 500000 > ${peaks_ps2}_masked_f &
	sort -k8,8nr ${peaks_ps3}_masked | head -n 500000 > ${peaks_ps3}_masked_f &
	sort -k8,8nr ${peaks_ps4}_masked | head -n 500000 > ${peaks_ps4}_masked_f &
	wait %1 %2 %3 %4
	sort -k8,8nr ${peaks_pool_ps1}_masked | head -n 500000 > ${peaks_pool_ps1}_masked_f &
	sort -k8,8nr ${peaks_pool_ps2}_masked | head -n 500000 > ${peaks_pool_ps2}_masked_f &
	wait %1 %2
	echo "done with management, time for idr"



	# IDR:
	# rep1 vs rep2, oracle: peaks from pooled sample
	${py3} ${idr} --input-file-type narrowPeak --rank p.value --soft-idr-threshold 0.1 -s ${peaks1}_masked_f ${peaks2}_masked_f -p ${peaks3}_masked_f --verbose -l ${current_wd}/comparison1.log -o ${current_wd}/${org}_${stage}_rep1_rep2.txt --plot &
	# rep1 pseudoreplicate1 Vs rep1 pseudoreplicate2, oracle: rep1 peaks
	${py3} ${idr} --input-file-type narrowPeak --rank p.value --soft-idr-threshold 0.1 -s ${peaks_ps1}_masked_f ${peaks_ps2}_masked_f -p ${peaks1}_masked_f --verbose -l ${current_wd}/comparison2.log -o ${current_wd}/${org}_${stage}_psr1_psr2.txt --plot &
	# rep2 pseudoreplicate1 Vs rep2 pseudoreplicate2, oracle: rep2 peakssudo surison3.log -o ${current_wd}/${org}_${stage}_psr3_psr4.txt --plot &
	# Pooled pseudoreplicates, oracle: peaks from pooled sample
	${py3} ${idr} --input-file-type narrowPeak --rank p.value --soft-idr-threshold 0.1 -s ${peaks_pool_ps1}_masked_f ${peaks_pool_ps2}_masked_f -p ${peaks3}_masked_f --verbose -l ${current_wd}/comparison4.log -o ${current_wd}/${org}_${stage}_poolPseudo.txt --plot &
	wait %1 %2 %3 %4
	echo 'DONE!'

	awk 'BEGIN{OFS="\t"} $12>='"1"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' ${current_wd}/${org}_${stage}_rep1_rep2.txt | sort | uniq | sort -k7n,7n  > ${TMPDIR}/results/${org}_${stage}_${samplesize}_idr01Peaks.bed

	# we need some cleanup:
	rm ${current_wd}/*.xls
	rm ${current_wd}/*masked*
	rm ${current_wd}/*bed
	rm ${current_wd}/*png

	tar vczf ${current_wd}.tar.gz ${current_wd}
	rm -r ${current_wd}

	turn=$((turn+1))

done

mkdir ${TMPDIR}/results/loops
mv ${TMPDIR}/results/*gz ${TMPDIR}/results/loops/

rsync -ah ${TMPDIR}/results/* ${from}/results_${org}_${stage}_${SLURM_JOB_ID}

rm -r $TMPDIR
