#!/home/ska/panos/venv/bin/python

from multiprocessing import Pool
import pysam
import sys
import pandas as pd
import gzip


def gen_chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        try:
            yield l[i:i+n]
        except StopIteration:
            yield l[i:]

def read_to_bedline(read, chrom):
    if read.is_reverse:
        s = "{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, max([read.pos -4, 0]), read.aend -4, ".",read.mapping_quality, "-")
    else:
        s = "{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, read.pos + 5 , read.aend + 5, ".",read.mapping_quality, "+")
    return s

def q_get_references(fpi):
    samfile = pysam.AlignmentFile(fpi, "rb")
    refs = samfile.references
    samfile.close()
    return refs

def one_worker_does(fpi, fpo,lor ):
    samfile = pysam.AlignmentFile(fpi, "rb")

    with gzip.open(fpo, "wb") as fo:
        for reference in lor:
            for read in samfile.fetch(reference):
                if read.is_proper_pair and read.tlen <= limit and read.mapping_quality >=10:
                    fo.write(read_to_bedline(read,reference))
    samfile.close()
    return
            

    
# The script expects 4 command line arguments:

# A file path to the input bam file
fpi = sys.argv[1]
# An upper limit for the fragment length.
limit = int(sys.argv[2])
# A number of processors for parallelization
p = int(sys.argv[3])
# A template for the resulting bed.gz files, something like : "./results/chunk_{}.bed.gz"
# this needs to include "{}" so that the script can name the various chunks
prefix = sys.argv[4]



# Some load balancing, we order the chromosomes by # of mapped reads and distribute
# chunks of chromosomes with relatively equal number of reads.
idx = pysam.idxstats(fpi)
idx = [x.rstrip().split("\t") for x in idx.split('\n')]

df = pd.DataFrame(idx, columns= ["chrom", "chrom_size","mapped", "unmapped"])
df = df[df['chrom']!='*']
df = df.iloc[:-1]
df.mapped = df.mapped.astype(int)
df = df.sort_values(by='mapped')

# careful, chrM gets ignored by default.
df = df[df['chrom']!='chrM']
sam = df.mapped.sum()
df['cum_mapped'] = df['mapped'].cumsum()

denom=sam/(int(p)*3)

df['yolo'] = (df.cum_mapped/(denom)).astype(int)

chrom_groups = df.groupby("yolo").apply(lambda g: g.chrom.tolist())

pool = Pool(processes=p)
for enum,chunk in enumerate(chrom_groups):
    fpo = prefix.format(enum)

    pool.apply_async(one_worker_does,args=(fpi, fpo, chunk,) )
    
pool.close()
pool.join()
del pool
