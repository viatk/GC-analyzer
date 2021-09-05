# GC Content Analyzer
Scripts that I wrote for analyzing GC content of genomic features. These scripts were extensively used during development of SS-seq and BP-seq methods published in [ Genome Research ](https://genome.cshlp.org/content/early/2021/06/24/gr.270082.120).

**Motivation:** Functional genomic elements such as, transcription start sites (TSS) or transcription end sites (TES), frequently have unusual GC content. It is important to take this into account when developing new NGS methods as many techniques themselves have AT/GC bias. The following scripts visualize GC content as a heatmap and as an average plot, estimate GC content bias of the functional elements of interest and generate random sequences of specified content for a control dataset.

**Usage:**

**To visualize GC content**

``Plot_GC_Content.R -- file_name worm_TAS.fasta --regionLength 2000 -- binSize 50 – numberOfHexbins 100 --label TES -- sortRegionSize 0 -- out_file_name worm_TAS ``

Below is GC content plotted for TESs in *C. elegans* using above command. As you can see from the plot, TESs are highly AT rich.
<img src=./TES_heatmap_plot.png> 
**To sort genomic features by GC content**

``Plot_GC_Content.R -- file_name worm_TAS.fasta --regionLength 2000 -- binSize 50 – numberOfHexbins 100 --label TSS -- sortRegionSize 400 -- out_file_name worm_TAS ``

Below is GC content plotted for TSSs in *C. elegans* using above command sorted by GC %. All TSSs have GC content higher than surrounding sequence.

<img src=./GC_Heatmap_Promoters_400.png> 

**To estimate GC bias**

``Plot_GC_Bias.R -- genome_file_name genomic_01.fa -- seq_file_name atac_coding_promoter.fasta –window 200``

Below is a plot comparing distribution of GC % of *C. elegans* TSSs to all the genomic sequences of the same length. *C. elegans TSSs* are significantly more GC rich than the rest of the genome.

<img src=./atac_coding_GC_bias.jpg> 

**To generate set of random sequences with GC% similar to genomic features**
