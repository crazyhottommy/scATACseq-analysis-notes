# scATAC-seq-analysis-notes
my notes for scATACseq analysis

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
* [Building gene regulatory networks from **single-cell** ATAC-seq and RNA-seq using Linked Self-Organizing Maps](https://www.biorxiv.org/content/early/2018/10/09/438937) github link https://github.com/csjansen/SOMatic
