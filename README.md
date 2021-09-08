# GC Content Analyzer
Scripts that I wrote for analyzing GC content of genomic features. These scripts were extensively used during development of SS-seq and BP-seq methods published in [ Genome Research ](https://genome.cshlp.org/content/early/2021/06/24/gr.270082.120).

**Motivation:** Functional genomic elements such as, transcription start sites (TSS) or transcription end sites (TES), frequently have unusual GC content. It is important to take this into account when developing new NGS methods as many techniques themselves have AT/GC bias. The following scripts visualize GC content as a heatmap and as an average plot, estimate GC content bias of the functional elements of interest and generate random sequences of specified content for a control dataset.

**Usage:**
Options:

        -f CHARACTER, --file=CHARACTER
                sequences file name

        -r INTEGER, --regionLength=INTEGER
                legnth of the region to be plotted

        -b INTEGER, --binSize=INTEGER
                bin size to be used in average plot

        -l CHARACTER, --label=CHARACTER
                label to put on the plots

        -s INTEGER, --sortRegionSize=INTEGER
                size of the region to use for sorting; if 0 sequences remain unsorted

        -n INTEGER, --numberOfHexbins=INTEGER
                number of hexbins to use in GC% heatmap; must be between 10 and 100

        -h, --help
                Show this help message and exit

**To visualize GC content**

``Rscript Plot_GC_Content.R --file worm_TAS.fasta --regionLength 2000 --binSize 50 --label "TSS" --sortRegionSize 0 --numberOfHexbins 100 ``

Below is GC content plotted for TESs in *C. elegans* using above command. As you can see from the plot, TESs are highly AT rich.
<img src=./TES_heatmap_plot.png> 
**To sort genomic features by GC content**

``Rscript Plot_GC_Content.R --file atac_coding_promoter.fasta --regionLength 2000 --binSize 50 --label "TSS" --sortRegionSize 300 --numberOfHexbins 100 ``

Below is GC content plotted for TSSs in *C. elegans* using above command sorted by GC %. All TSSs have GC content higher than surrounding sequence.

<img src=./GC_Heatmap_Promoters_400.png> 

**To estimate GC bias**

``Plot_GC_Bias.R -- genome_file_name genomic_01.fa -- seq_file_name atac_coding_promoter.fasta –window 200``

Below is a plot comparing distribution of GC % of *C. elegans* TSSs to all the genomic sequences of the same length. *C. elegans TSSs* are significantly more GC rich than the rest of the genome.

<img src=./atac_coding_GC_bias.jpg> 

**To generate set of random sequences with GC% similar to genomic features**
