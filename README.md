# Mapping ATACseq reads
We receive the paired-end reads from the sequencing facility with their adaptors already trimmed. Let the two raw fastq files be reads_1.fq and reads_2.fq . We map the reads on the reference genome and end up with a sorted and indexed bam file with the following commands:

```sh
>> bowtie2 -x ${index} -1 reads_1.fq -2 reads_2.fq --very-sensitive -X 2000 -I 0  | \   # Map with bowtie
samtools view -Shu - | \                                                                # pipe sam to bam
samtools sort -T ${TMPDIR}/temp - | \                                                   # sort the bam file 
samtools rmdup -s - - > ${rbam}                                                         # remove duplicates

>> samtools index ${rbam}                                                    # index the resulting bam file
```


# Filtering the nucleosome-free reads, +4-5 ATACseq correction
We want to only use reads whose mate is less than X bp (in this work X=120) away. This ensures that we deplete our
input of transposase events that happened on both sides of a nucleosome, and thus allows us to get higher confidense
peaks. 

We use a home-made python script which uses the pysam library to read the bam file containing all mapped reads.
The script keeps only reads considered to be "proper pair" by the mapper (equivalent to samtools view -f 2),
with a mapping quality >=10 (sam file column #5) and with |fragment_length| <= 120 ( sam file column #9).

This script can be found in this repository as *filter_nf_reads.py*
