    
# Comparing WGCNA clusters from the two species    

### The clusters
We use WGCNA to cluster genes by their RNA activity in various TRNAseq experiments, both in developmental stages and pure tissue samples.  This analysis assigns each gene in a single cluster.  These clusters beg names like "blue". "green" etc.    
We  manually explore the clusters, and by determining at which experiment the genes of the cluster seem to be differentially expressed (but also by GO analysis on each cluster), we assign more meaningful names to the clusters , like "liver", "muscle" etc.

### Other data
* ##### Basal Regions 
    We take each gene's TSS, and assign from 5kb upstram to 1kb downstream as the gene's BASAL region.
* ##### Peaks
    Peaks have been called with the IDR analysis in each individual ATAC sample. For this analysis we use a merged set of peaks form all ATAC experiments per organism. That means that each gene will be associated with all peaks that we can possibly associate it with from our ATACs.
* ##### PWMs
    First we map the vertebrate.amphioxus.clustered.v2.00.filtered.pwm, Then they are merged even more, into "super families"
    so in the end each "TF" considered for enrichment is in fact a superfamily/superset of clustered pwms.    
     
### The analysis
1. We associate WGCNAclusters --> genes --> BASAL regions --> ATAC peaks --> PWM hits
2. The core for this analysis then, is a table that looks like this:
    |index|peakID|geneID|cluster|pwm|
    |---|---|---|---|---|
    |1|0|BL09450|blue|851;876;1065|
    |5|0|BL09450|blue|66;15742;15743;18247|
    |6|0|BL09450|blue|66;15742;15743;18247|
    |7|0|BL09450|blue|15284;15676;21479|
    |8|0|BL09450|blue|583|
3.  We can now ask the following questions (I can do a 'gene-based enrichment, peak-based or pwm-based) to calculate an enrichment per cluster per PWM.
    3.1 Population Size (**M**)
    *   total \# of genes 
    *   **or**  total \# of peaks
    *   **or**  total \# of pwm hits

    3.2 Success in population (for each PWM (or PWM superset)) (**N**)
    *   How many genes have at least one hit of this pwm?
    *   **or** how many peaks have a hit for this PWM? (in all genes)
    *   **or** how many hits of the PWM did I get in total (in all peaks, in all genes)
    
    3.3 Sample size: (**n**)
    * How many genes in this cluster?
    * How many ATACpeaks in this cluster?
    * How many PWM hits in this cluster?
    
    3.4 Success in sample (**x**)
    * How many genes in this cluster have a hit of this PWM?
    * How many peaks in this cluster have a hit of this PWM?
    * How many hits of this PWM do I have in this cluster?

4. Based on the above numbers, and with scipy's [hypergeometric "survival function"](https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.hypergeom.html)
    We calculate an "enrichment pvalue" per cluster-PWM pair ( **hypergeom.sf(x, M, n, N)**)
    Or alternatively **enrichment = (x/M) / ( (N/M) * (n/M) )** (as given to me by george )

5. Now we can calculate the correlation between two clusters, based on their enrichment values for all PWMs
6. Finally we throw this matrix of correlations into a clustermap function that gives us the final heatmap
    of *"cluster similarity based on their regulatory content"*


    
    
    































