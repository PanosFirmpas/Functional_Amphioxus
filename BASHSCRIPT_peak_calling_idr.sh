

###########################################################################################################################################################################################################################
# SOME NOTES:
#
# This script is for illustration purposes only. It should run well with few modifications, mostly having to do with paths
# DONT EXPECT THIS TO RUN WITHOUT ANY CHANGES
# 
# If you need help runing this, please feel free to contact me : panosfirbas [at] gmail [dot] com
# or in the issues section of this github repository
#
# Parallelization:
# we send the command to the background with & so we can call more commands in the same bash script,
# then wait for the commands to finish with 'wait'
# We parallelize up to 4 commands in this script
#
###########################################################################################################################################################################################################################


###########################################################################################################################################################################################################################
# SOFTWARE USED IN THIS SCRIPT
#
# If you have the following already in your $PATH, then you can replace all ${macs} instances with macs and ${idr} with idr
# we don't so we need to point to them.
macs2="/path/to/bin/macs2"	# the macs2 peak caller: https://github.com/taoliu/MACS
#
#
# This is the python with which idr has been installed, macs2 requires python2 so we have to get a bit creative and explicitely tell the system which python to use for idr (py2 is default)
py3=/path/to/executable/of/python3   
idr='/path/to/idr'		         # the python idr package: https://github.com/nboley/idr
#
# bedtools : http://bedtools.readthedocs.io
# zcat , shuf, awk, cat, sort, uniq, wc, head, tail etc are all standard GNU/Linux Command-Line Tools/features
# If you are not in such a system, please consider switching to non-proprietory Free and open-source software.
###########################################################################################################################################################################################################################


###########################################################################################################################################################################################################################
# COMMAND LINE ARGUMENTS
#
# Organism and stage NAMES, just so the script is not hardcoded, they have no computational significance
org=$1			# for example: 'zebra'
stage=$2		# for example: 'dome'
# the following are important:
gsize=$3		# The genome size of the organism where we're peak calling, very important for macs2 background calculation
mask=$4			# A bed file with regions of the genome that are to be ignored, repeated elements etc. Taken from ucsc (?)
###########################################################################################################################################################################################################################



###########################################################################################################################################################################################################################
## Some pre-processing before we peak-call
###########################################################################################################################################################################################################################
TMPDIR='/path/to/where/you/will/work'    # we run this on a cluster with a job scheduler so this is set to /scratch/temp_${SLURM_JOB_ID}/
mkdir ${TMPDIR}/results
mkdir ${TMPDIR}/data
cd ${TMPDIR}

# The nucleosome free reads that have been precomputed
# the filter_nf_reads.py script gives a number of bed files per experiment:
# replicate_1_nf_reads_1.gz, replicate_1_nf_reads_2.gz .... replicate_1_nf_reads_N.gz
replic1=/path/to/${org}_${stage}_replicate_1_nf_reads_*.gz
replic2=/path/to/${org}_${stage}_replicate_2_nf_reads_*.gz

# A little cleaning of the NF reads,  no Start <=0, no ChrM, no low quality, no more than 60% overlap with masked elements
zcat ${replic1} | awk '( $2>=0 ) && ( $5 > 10 ) && ( $1 != chrM )' | bedtools intersect -v -a - -b ${mask} -f 0.6 | shuf > ./data/replic1_clean.bed &
zcat ${replic2} | awk '( $2>=0 ) && ( $5 > 10 ) && ( $1 != chrM )' | bedtools intersect -v -a - -b ${mask} -f 0.6 | shuf > ./data/replic2_clean.bed &
wait %1 %2

# count reads in the two replicates 
replic1_lines=$( cat ./data/replic1_clean.bed | wc -l )
replic2_lines=$( cat ./data/replic2_clean.bed | wc -l )
# we'll get half for each pseudoreplicate
replic1_sample_size=$(( (replic1_lines + 1) / 2 ))
replic2_sample_size=$(( (replic2_lines + 1) / 2 ))
#Make the replicate 1 pseudoreplicates
head -n ${replic1_sample_size} ./data/replic1_clean.bed > ./data/replic1_psr1.bed
tail -n ${replic1_sample_size} ./data/replic1_clean.bed > ./data/replic1_psr2.bed
#Make the replicate 2 pseudoreplicates
head -n ${replic2_sample_size} ./data/replic2_clean.bed > ./data/replic2_psr1.bed
tail -n ${replic2_sample_size} ./data/replic2_clean.bed > ./data/replic2_psr2.bed

# Then, the Pooled pseudoreplicates
pooled_lines=$(( replic1_lines + replic2_lines ))
pool_sample_size=$(( (pooled_lines + 1) / 2 ))
cat ./data/replic1_clean.bed ./data/replic2_clean.bed | shuf > ./data/pooled.bed
head -n ${pool_sample_size} ./data/pooled.bed > ./data/pool_psr1.bed
tail -n ${pool_sample_size} ./data/pooled.bed > ./data/pool_psr2.bed

###########################################################################################################################################################################################################################
## THE PEAK CALLING!
###########################################################################################################################################################################################################################

# replicate1, replicate 2, and pooled reads
${macs2} callpeak --nomodel --keep-dup 1 --llocal 10000 --extsize 74 --shift -37  -p 0.07 -f BED -g ${gsize} -t ./data/replic1_clean.bed -n ${org}_${stage}_rep1_nf --outdir ${TMPDIR}/results/ &
${macs2} callpeak --nomodel --keep-dup 1 --llocal 10000 --extsize 74 --shift -37  -p 0.07 -f BED -g ${gsize} -t ./data/replic2_clean.bed -n ${org}_${stage}_rep2_nf --outdir ${TMPDIR}/results/ &
${macs2} callpeak --nomodel --keep-dup 1 --llocal 10000 --extsize 74 --shift -37  -p 0.07 -f BED -g ${gsize} -t ./data/replic1_clean.bed ./data/replic2_clean.bed -n ${org}_${stage}_comb_nf --outdir ${TMPDIR}/results/ &
wait %1 %2 %3

#Peak call for the replicate 1 pseudoreplicates
${macs2} callpeak --nomodel --keep-dup 1 --llocal 10000 --extsize 74 --shift -37  -p 0.07 -f BED -g ${gsize} -t ./data/replic1_psr1.bed -n ${org}_${stage}_rep1_psr1 --outdir ${TMPDIR}/results/ &
${macs2} callpeak --nomodel --keep-dup 1 --llocal 10000 --extsize 74 --shift -37  -p 0.07 -f BED -g ${gsize} -t ./data/replic1_psr2.bed -n ${org}_${stage}_rep1_psr2 --outdir ${TMPDIR}/results/ &
#Peak call for the replicate 2 pseudoreplicates
${macs2} callpeak --nomodel --keep-dup 1 --llocal 10000 --extsize 74 --shift -37  -p 0.07 -f BED -g ${gsize} -t ./data/replic2_psr1.bed -n ${org}_${stage}_rep2_psr1 --outdir ${TMPDIR}/results/ &
${macs2} callpeak --nomodel --keep-dup 1 --llocal 10000 --extsize 74 --shift -37  -p 0.07 -f BED -g ${gsize} -t ./data/replic2_psr2.bed -n ${org}_${stage}_rep2_psr2 --outdir ${TMPDIR}/results/ &
wait %1 %2 %3 %4

#Peak calling for the pooled pseudoreplicates
${macs2} callpeak --nomodel --keep-dup 1 --llocal 10000 --extsize 74 --shift -37  -p 0.07 -f BED -g ${gsize} -t ./data/pool_psr1.bed -n ${org}_${stage}_pool_psr1 --outdir ${TMPDIR}/results/ &
${macs2} callpeak --nomodel --keep-dup 1 --llocal 10000 --extsize 74 --shift -37  -p 0.07 -f BED -g ${gsize} -t ./data/pool_psr2.bed -n ${org}_${stage}_pool_psr2 --outdir ${TMPDIR}/results/ &
wait %1 %2

###########################################################################################################################################################################################################################
## Some pre-processing before IDR
###########################################################################################################################################################################################################################

# define some names for ease
peaks1=./results/${org}_${stage}_rep1_nf_peaks.narrowPeak
peaks2=./results/${org}_${stage}_rep2_nf_peaks.narrowPeak
peakspool=./results/${org}_${stage}_comb_nf_peaks.narrowPeak
peaks_ps11=./results/${org}_${stage}_rep1_psr1_peaks.narrowPeak
peaks_ps12=./results/${org}_${stage}_rep1_psr2_peaks.narrowPeak
peaks_ps21=./results/${org}_${stage}_rep2_psr1_peaks.narrowPeak
peaks_ps22=./results/${org}_${stage}_rep2_psr2_peaks.narrowPeak
peaks_pool_ps1=./results/${org}_${stage}_pool_psr1_peaks.narrowPeak
peaks_pool_ps2=./results/${org}_${stage}_pool_psr2_peaks.narrowPeak

# Clean with rep mask. This time we are removing entire peaks so we set a less stringent threshold of 0.33
# if more than 0.33 of a peak overlaps a masked element, the peak is discarded
bedtools intersect -v -a ${peaks1} -b ${mask} -f 0.33 > ${peaks1}_masked &
bedtools intersect -v -a ${peaks2} -b ${mask} -f 0.33 > ${peaks2}_masked &
bedtools intersect -v -a ${peakspool} -b ${mask} -f 0.33 > ${peakspool}_masked &
wait %1 %2 %3
bedtools intersect -v -a ${peaks_ps11} -b ${mask} -f 0.33 > ${peaks_ps11}_masked &
bedtools intersect -v -a ${peaks_ps12} -b ${mask} -f 0.33 > ${peaks_ps12}_masked &
bedtools intersect -v -a ${peaks_ps21} -b ${mask} -f 0.33 > ${peaks_ps21}_masked &
bedtools intersect -v -a ${peaks_ps22} -b ${mask} -f 0.33 > ${peaks_ps22}_masked &
wait %1 %2 %3 %4
bedtools intersect -v -a ${peaks_pool_ps1} -b ${mask} -f 0.33 > ${peaks_pool_ps1}_masked &
bedtools intersect -v -a ${peaks_pool_ps2} -b ${mask} -f 0.33 > ${peaks_pool_ps2}_masked &
wait %1 %2


###########################################################################################################################################################################################################################
## THE IDR !!
###########################################################################################################################################################################################################################
# rep1 vs rep2, oracle: peaks from pooled sample
${py3} ${idr} -i 0.01 -s ${peaks1}_masked ${peaks2}_masked -p ${peakspool}_masked --verbose -l ./results/comparison1.log -o ./results/${org}_${stage}_rep1_rep2.txt --plot &
# rep1 pseudoreplicate1 Vs rep1 pseudoreplicate2, oracle: rep1 peaks
${py3} ${idr} -i 0.01 -s ${peaks_ps11}_masked ${peaks_ps12}_masked -p ${peaks1}_masked --verbose -l ./results/comparison2.log -o ./results/${org}_${stage}_psr1_psr2.txt --plot &
# rep2 pseudoreplicate1 Vs rep2 pseudoreplicate2, oracle: rep2 peaks
${py3} ${idr} -i 0.01 -s ${peaks_ps21}_masked ${peaks_ps22}_masked -p ${peaks2}_masked --verbose -l ./results/comparison3.log -o ./results/${org}_${stage}_psr3_psr4.txt --plot &
# Pooled pseudoreplicates, oracle: peaks from pooled sample
${py3} ${idr} -i 0.01 -s ${peaks_pool_ps1}_masked ${peaks_pool_ps2}_masked -p ${peakspool}_masked --verbose -l ./results/comparison4.log -o ./results/${org}_${stage}_poolPseudo.txt --plot &
wait %1 %2 %3 %4

# Filter the peaks with idr <=0.01. This is our proper peak list, the other idrs were for diagnostic reasons.
awk 'BEGIN{OFS="\t"} $12>='"1"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' ./results/${org}_${stage}_rep1_rep2.txt | sort | uniq | sort -k7n,7n  > ./results/${org}_${stage}_idr01Peaks.bed


###########################################################################################################################################################################################################################
# Final note: You can now discard the ${TMPDIR}/data folder and keep the ${TMPDIR}/results/
###########################################################################################################################################################################################################################
