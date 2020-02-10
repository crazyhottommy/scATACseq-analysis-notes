# scATAC-seq-analysis-notes
my notes for scATACseq analysis


### paper to read

[Prospective, brain-wide labeling of neuronal subclasses with enhancer-driven AAVs](https://www.biorxiv.org/content/early/2019/01/20/525014?rss=1) Neuroscientists need ever-more specific methods to restrict expression of tools to defined neuronal subpopulations. This amazing study from @AllenInstitute uses scATAC-seq & scRNA-seq to discover enhancers that are even more specific than driver lines.

### ATAC-seq QC

### protocols 

[An improved ATAC-seq protocol reduces background and enables interrogation of frozen tissues](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5623106/)

#### Fragment length distribution

A blog post by [Xi Chen](https://dbrg77.wordpress.com/2017/02/10/atac-seq-insert-size-plotting/)

>The successful construction of a ATAC library requires a proper pair of Tn5 transposase cutting events at the ends of DNA. In the nucleosome-free open chromatin regions, many molecules of Tn5 can kick in and chop the DNA into small pieces; around nucleosome-occupied regions, Tn5 can only access the linker regions. Therefore, in a normal ATAC-seq library, you should expect to see a sharp peak at the <100 bp region (open chromatin), and a peak at ~200bp region (mono-nucleosome), and other larger peaks (multi-nucleosomes)

https://github.com/GreenleafLab/NucleoATAC/issues/18

>there might be some artifact with how the aligner deals with fragments where the forward read and reverse read are exact reverse complements of each other. I know that Bowtie (1 but not 2) has some issue with those reads.

### ATAC-seq

>Some may notice that the peaks produced look both like peaks produced from the TF ChIP-seq pipeline as well as the histone ChIP-seq pipeline. This is intentional, as ATAC-seq data looks both like TF data (narrow peaks of signal) as well as histone data (broader regions of openness).

* [ATACseqQC ](http://bioconductor.org/packages/release/bioc/html/ATACseqQC.html) a bioconductor package for quality control of ATAC-seq data.
* [RASQUAL](https://github.com/dg13/rasqual) (Robust Allele Specific QUAntification and quality controL) maps QTLs for sequenced based cellular traits by combining population and allele-specific signals. [paper: Fine-mapping cellular QTLs with RASQUAL and ATAC-seq](http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3467.html) 
* [ATAC-seq Forum](https://sites.google.com/site/atacseqpublic/home?pli=1)  
* [Single-cell ATAC-Seq](http://cole-trapnell-lab.github.io/projects/sc-atac/) 
* [A rapid and robust method for single cell chromatin accessibility profiling](https://www.biorxiv.org/content/early/2018/04/27/309831)
* [Global Prediction of Chromatin Accessibility Using RNA-seq from Small Number of Cells](http://www.biorxiv.org/content/early/2016/01/01/035816)  from RNA-seq to DNA accessibility. [tool on github](https://github.com/WeiqiangZhou/BIRD)  
* [NucleoATAC](https://github.com/GreenleafLab/NucleoATAC):Python package for calling nucleosomes using ATAC-Seq data 
* [chromVAR: Inferring transcription factor variation from **single-cell** epigenomic data](http://biorxiv.org/content/early/2017/02/21/110346) scATAC-seq
* [ENCODE ATAC-seq guidelines](https://www.encodeproject.org/data-standards/atac-seq/)
* [Brockman](https://carldeboer.github.io/brockman.html) is a suite of command line tools and R functions to convert genomics data into DNA k-mer words representing the regions associated with a chromatin mark, and then analyzing these k-mer sets to see how samples differ from each other. This approach is primarily intended for **single cell** genomics data, and was tested most extensively on single cell ATAC-seq data
* [Reproducible inference of transcription factor footprints in ATAC-seq and DNase-seq datasets via protocol-specific bias modeling](https://www.biorxiv.org/content/early/2018/03/19/284364)
* [msCentipede](http://rajanil.github.io/msCentipede/) is an algorithm for accurately inferring transcription factor binding sites using chromatin accessibility data (Dnase-seq, ATAC-seq) and is written in Python2.x and Cython.
* The Differential ATAC-seq Toolkit [(DAStk)](https://biof-git.colorado.edu/dowelllab/DAStk) is a set of scripts to aid analyzing differential ATAC-Seq data.
* [Identification of Transcription Factor Binding Sites using ATAC-seq](https://www.biorxiv.org/content/early/2018/07/17/362863)  We propose HINT-ATAC, a footprinting method that addresses ATAC- seq specific protocol artifacts
* [HMMRATAC](https://github.com/LiuLabUB/HMMRATAC) splits a **single** ATAC-seq dataset into nucleosome-free and nucleosome-enriched signals, learns the unique chromatin structure around accessible regions, and then predicts accessible regions across the entire genome. We show that HMMRATAC outperforms the popular peak-calling algorithms on published human and mouse ATAC-seq datasets.
* [ChromA](https://github.com/marianogabitto/ChromA): Chromatin Landscape Annotation Tool. 
ChromA is a probabilistic model to annotate chromatin regions into accessible or inaccessible, open or closed, based on their ATACseq profile. ChromA can process bulk datasets, single-cell or integrate information from a combination of both. Even more, ChromA can integrate information from different replicates or different cellular populations to create a consensus representation of chromatin accessibility.

### peak calling

* [MACS2](https://github.com/taoliu/MACS/issues/145)

> --shift -100 --extsize 200 will amplify the 'cutting sites' enrichment from ATAC-seq data. So in the end, the 'peak' is where Tn5 transposase likes to attack. The fact is that, although many information such as the insertion length and the other mate alignment is ignored, such result is still usable. Especially when the short fragment population is extremely dominant, the final output won't be off much.

```bash
macs2 --nomodel --keepdup all --shift -100 --extsize 200
macs2 -f BAMPE
```
* [generich](https://informatics.fas.harvard.edu/atac-seq-guidelines.html#peak) written by John, previous labmates at Harvard FAS informatics. will take a look!

### motif analysis

* [marge](https://marge.netlify.com/)An API for Analysis of Motifs Using HOMER in R:  https://www.biorxiv.org/content/10.1101/249268v1
* [homerkit](https://github.com/slowkow/homerkit)

### Dimension Reduction

* [clustering scATACseq data: the TF-IDF way](https://divingintogeneticsandgenomics.rbind.io/post/clustering-scatacseq-data-the-tf-idf-way/) My effort to replicate some of the studies.
* [Dimensionality Reduction for scATAC Data](http://andrewjohnhill.com/blog/2019/05/06/dimensionality-reduction-for-scatac-data/) A blog post by Andrew Hill. Very informative.

### clustering

* [cisTopic: cis-regulatory topic modeling on single-cell ATAC-seq data](https://www.nature.com/articles/s41592-019-0367-1) github https://github.com/aertslab/cistopic

### footprint
* [HINT](https://www.regulatory-genomics.org/) tutorial: https://www.regulatory-genomics.org/hint/tutorial/
paper: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1642-2 single cell ATAC tutorial is on the way.
* [seqoutBias](https://guertinlab.github.io/seqOutBias_Vignette/part3.html) remove Tn5 cut bias.
* [scOpen: chromatin-accessibility estimation of single-cell ATAC data](https://www.biorxiv.org/content/biorxiv/early/2019/12/05/865931.full.pdf) tutorial https://www.regulatory-genomics.org/hint/tutorial-differential-footprints-on-scatac-seq/

### piplines

* [scATAC-pro: a comprehensive workbench for single-cell chromatin accessibility sequencing data](https://www.biorxiv.org/content/10.1101/824326v2) https://github.com/wbaopaul/scATAC-pro



### integrate scATAC and scRNAseq

* [MOFA+: analysis of matching scRNA-seq and scATAC-seq data](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/SNARE_seq.html)
* Buenrostro et al., Integrated Single-Cell Analysis Maps the Continuous Regulatory Landscape of Human Hematopoietic Differentiation, Cell (2018), https://doi.org/10.1016/j.cell.2018.03.074
* Duren et al., Integrative analysis of single-cell genomics data by coupled nonnegative matrix factorizations, PNAS (2018), https://doi.org/10.1073/pnas.1805681115
* [Prospective, brain-wide labeling of neuronal subclasses with enhancer-driven AAVs](https://www.biorxiv.org/content/10.1101/525014v2)
* [Integrative single-cell analysis of transcriptional and epigenetic states in the human adult brain](https://www.nature.com/articles/nbt.4038) Jean Fan kindly shared her scripts for the paper at https://github.com/JEFworks/Supplementary-Code/tree/master/snDropSeq_scTHSseq. will need to ask her when I have questions.
* [Building gene regulatory networks from **single-cell** ATAC-seq and RNA-seq using Linked Self-Organizing Maps](https://www.biorxiv.org/content/early/2018/10/09/438937) github link https://github.com/csjansen/SOMatic
* [Comprehensive integration of single cell data](https://www.biorxiv.org/content/early/2018/11/02/460147) integrating scATAC and scRNAseq in the [`seurat`](https://www.satijalab.org/seurat) v3 R package.
* [Integrative analysis of single cell genomics data by coupled nonnegative matrix factorizations](http://web.stanford.edu/~zduren/CoupledNMF/)
* Multi-Omics Factor Analysis (MOFA) http://bioconductor.org/packages/release/bioc/html/MOFA.html integrate scATAC and scRNAseq?
* [Garnett](https://cole-trapnell-lab.github.io/garnett/) automatic annotation of scRNAseq and scATACseq data sets.
* [Fluent genomics with plyranges and tximeta](https://sa-lee.github.io/fluentGenomics/articles/fluentGenomics.html) integrate bulk RNAseq and ATACseq data using bioconductor packages and using plyranges (dplyr for GRanges)!

### predicting ATAC peak target gene

* [Cicero](https://cole-trapnell-lab.github.io/cicero-release/)
* [ChIA-Drop: Multiplex chromatin interactions with single-molecule precision](https://www.nature.com/articles/s41586-019-0949-1) by Yijun Ruan at JAX lab. His office is next to my previous postdoc advisor Roel Verhaak's.

### pipelines

* [SnapATAC](https://github.com/r3fang/SnapATAC) much faster than `cellranger-atac`.
