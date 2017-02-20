# 1. Mapping ATACseq reads
We receive the paired-end reads from the sequencing facility with their adaptors already trimmed. Let the two raw fastq files be reads_1.fq and reads_2.fq (for paired-end sequencing). We map the reads on the reference genome and end up with a sorted and indexed bam file with the following commands:

```sh
>> bowtie2 -x ${bt2_index} -1 reads_1.fq -2 reads_2.fq --very-sensitive -X 2000 -I 0  | \   # Map with bowtie
    samtools view -Shu - | \                                                                # pipe sam to bam
    samtools sort - | \                                                                     # sort the bam file
    samtools rmdup -s - - > ${rbam}                                                         # remove duplicates

>> samtools index ${rbam}                                                    # index the resulting bam file
```

**${bt2_index}** : the appropriate bowtie2 index name (depends on how you have "installed" the organism's genome in bowtie2)
**${rbam}**  : the name of the resulting bam file, something like amphi_15hpf_rep1.bam


# 2. Getting the nucleosome-free reads, +5-4 ATACseq correction
We want to only use reads whose mate is less than X bp (in this work X=120) away. This ensures that we deplete our
input of transposase events that happened on both sides of a nucleosome, and thus allows us to get higher confidense
peaks. 

We use a home-made python script which uses the pysam library to read the bam file containing all mapped reads.
The script keeps only reads considered to be "proper pair" by the mapper (equivalent to samtools view -f 2),
with a mapping quality >=10 (sam file column #5) and with |fragment_length| <= 120 ( sam file column #9).
Reads mapped to the + strand are moved by +5 bp and reads mapping to the - strand by -4bp.

The script is parallelized and outputs a number of gzipped bed files where each bedline is one read, with
its mapping quality in the score column.

This script can be found in this repository as *00.2_filter_nf_reads.py*

```sh
# Example usage:
>> python 00.2_filter_nf_reads.py ${rbam} ${fs_limit} ${processors} ${prefix}
```
**${fs_limit}** : The upper FragmentSize limit, only reads with a fragment size <= this will be output.
**${processors}** : How many processors you want to employ 
**${prefix}**   : A string which contains "{}", for example: "amphi_15hpf_rep1_nfreads_{}.bed.gz". The script
will output a gzipped file per processor used and use this prefix to name the output files by replacing the {}
with a unique number. You can concatenate the resulting files with zcat 
(i.e. 
```
>> zcat amphi_15hpf_rep1_nfreads_*.bed.gz > amphi_15hpf_rep1_nfreads.bed.gz 
```
)





# 3. Peak calling
Please see 00.3_peak_calling_idr.sh for a detailed script of peak calling and idr. 
*Briefly*, we used [macs2](https://github.com/taoliu/MACS) for low-thresholded peak-calling 
```sh
>> macs2 callpeak --nomodel --keep-dup 1 --llocal 10000 --extsize 74 --shift -37  -p 0.07 -g ${gsize} \
    -t our_nucFree_reads_forRep1_fromStep2.bed -n repl_1_peaks.bed
```
Then, we use the [IDR framework](https://github.com/nboley/idr) to take advantage of our replicates.
>In layman's terms, the IDR method compares a pair of ranked lists of identifications (such as ChIP-seq peaks). These ranked lists should not be pre-thresholded i.e. they should provide identifications across the entire spectrum of high confidence/enrichment (signal) and low confidence/enrichment (noise). The IDR method then fits the bivariate rank distributions over the replicates in order to separate signal from noise based on a defined confidence of rank consistency and reproducibility of identifications i.e the IDR threshold.    
    

```sh
# Please see 00.3_peak_calling_idr.sh for a detailed script of peak calling and idr. 
>> idr -i 0.1 -s repl_1_peaks.bed repl_2_peaks.bed -p pooled_repl_peaks.bed -o idr_out.txt
>> awk 'BEGIN{OFS="\t"} $12>='"1"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' idr_out.txt uniq  | \
    sort -k7n,7n  > idr01Peaks.bed
```


# 4. Mapping PWMs to ATACseq peaks
We usea large set of clustered vertebrate and amphioxus PWMs and the [gimmemotifs](http://gimmemotifs.readthedocs.io/) tool to determine thresholds per PWM and then to map them on the genome.  
    
* #####  Determining a threshold for each motif
    Different PWMs have different maximum and minimum possible scores and thus different propensities to
    give a 'hit' on any given sequence. We want to be stricter with small/promiscuous PWMs and more lenient 
    with larger ones.

    'gimme threshold' calculates an optimal threshold per PWM, based on their propensity to bind on a given control
    sequence. We will generate a fake background sequence with 'gimme background', based on the gc content of non-ATACseq
    regions of our genome and use it to calculate the PWM thresholds.
    
    1. Genomic background:  
    Let 'peaks.bed' be the ATACseq peaks as determined by our previous analysis.    
        
        ```
        >> bedtools shuffle -i peaks.bed -g genome_of_interest.txt > real_background.bed  
        >> bedtools getfasta -fi danRer10.fa -bed real_background.bed > real_background.fa  
        >> gimme background -i real_background.fa -l 500 -n 50000 -f FASTA gc_background.fa gc  
        ```  

    2. Gimme threshold  
        Let 'factors.pwm' be the txt file that contains our pwms of interest  
        
        ```  
        >> gimme threshold factors.pwm gc_background.fa 0.01 > factors_thresholds.txt  
        ```  
        will give us thresholds with 0.01fdr stringency.  

* #####  Mapping the pwms on ATACseq peaks  
    We can now map the PWMs with their appropriate threshold on the real ATACseq regions:  
    ```  
    >>> bedtools getfasta -fi danRer10.fa -bed peaks.bed > peaks.fa
    >>> gimme scan -b -c factors_thresholds.txt peaks.fa factors.pwm > factors_mapped.bed
    ```






