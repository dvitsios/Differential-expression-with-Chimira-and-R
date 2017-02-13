# Chimira website:
<http://wwwdev.ebi.ac.uk/enright-dev/chimira>

## Tutorial #
*  You can use Chimira (Other tools tab) to infer the 3’ adapter from your FASTQ files:
<br/>
<http://wwwdev.ebi.ac.uk/enright-dev/chimira>

*  Then, you need to use the Clean & Run tab of Chimira in order to clean your files from the 3’ adapter and retrieve the expression counts.
In case you have e.g. 2 control and 2 treatment samples you can upload them altogether in order to get their counts in a single file in the end.

*  As soon as the analysis is complete, you can check the quality of the trimming from the ‘QC plots’ tab on the results page.
The most important thing there is to make sure that in the ‘Read lengths distribution after trimming’ plot the majority of your reads are around 21-23nt long.
This means that your cleaned samples will now contain mainly miRNA reads.

*  The only other result that you need to use for the rest of the analysis is the “Plain counts per file (Raw)” file, which can be found on a panel at the top-right corner of the results page.


*  After downloading this file, you just need to rename the first line so that it reflects the samples that you are analysing. For instance, you can rename the default line:

```
files   1/  2/  3/  4/

into

files   CTL1    CTL2    KO1 KO2
```

Please bear in mind that the required formatting for the renamed columns is :

nameA1  nameA2  nameB1  nameB2 &nbsp;&nbsp;&nbsp;<i>(tab separated)</i>


where ‘nameA' and ‘nameB’ can be any names descriptive of your samples and the suffixes ‘1’ and ‘2’ just denote the two replicates used for the analysis.

*  Now you are ready to analyse the file that you have downloaded (let’s rename it into ‘counts_clean.txt’) using the R script (analyse_counts.R):

The script accepts 3 input arguments: 
* the counts file
* the descriptive name of the first sample type
* the descriptive name of the second sample type

You can call the R scirpt from console e.g. like that:

```

Rscript analyse_counts.R counts_clean.txt CTL KO

```

..* reminder: in that case your counts file should be named ‘counts_clean.txt’ and the first line of that file should be:
<br/>
files   CTL1    CTL2    KO1 KO2



## Results ##

* This script is normalising your count data with [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and plots some basic QC at first (file Basic_QC.pdf), which contains a heatmap with sample-to-sample correlation, PCA analysis, dispersion plot and barplots with counts before and after normalisation.

* In the ‘Differential_Expression_analysis.pdf’ plot you get the differential expression scatterplot between the samples (after merging the replicates from each sample).
In this plot, all the dots highlighted in light-blue demonstrate statistically significant differential expression ( absolute(fold_change) >  2 & p-value < 0.05 ).
Among those, the let7-family miRNAs are also highlighted on the same plot.

* Apart from that, in the same pdf file, you get a scatterplot of diff expression of let7 miRNAs between your samples as well as barplots with the let7-family expression across all replicates and/or merged samples.

* Finally, some extra .txt files are being created which contain your normalised count data, a list of differentially expressed miRNAs (with the log2 fold change and adjusted p-value) and also the let7 counts (normal and log2 scale).

<br/>
An example with an input file and the results can be found in the R_diff_expr_analysis-example/ directory 
