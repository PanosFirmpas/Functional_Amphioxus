# Mapping ATACseq reads
We receive the paired-end reads from the sequencing facility with their adaptors already trimmed. Let the two raw fastq files be reads_1.fq and reads_2.fq . We map the reads on the reference genome and end up with a sorted and indexed bam file with the following commands:

```sh
>> bowtie2 -x ${index} -1 reads_1.fq -2 reads_2.fq --very-sensitive -X 2000 -I 0  | \   # Map with bowtie
samtools view -Shu - | \                                                                # pipe sam to bam
samtools sort -T ${TMPDIR}/temp - | \                                                   # sort the bam file 
samtools rmdup -s - - > ${rbam}                                                         # remove duplicates

>> samtools index ${rbam}                                                       # index the resulting bam file
```





