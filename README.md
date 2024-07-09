# Yeast selective mitoribosome profiling analysis 

This repository includes the custom scripts, some downloaded programs, and many sequence files needed to analyze yeast mitoribosome profiling data and make publication-quality plots. Note the reference genome will have to be downloaded separately. The commands as written are to be used with the SLURM job scheduler and many paths throughout will need to be updated for your directory structure. Additionally several packages/modules need to be loaded or in your path; e.g. python, samtools, bedtools, cutadapt, R, bowtie, bamtools, fastqc, java, star, picard .

Listed below is an overview of the functions provided. Examples for the commands to be used to accomplish each step are provided in 'yMitoRP_pipelineCommands_README.txt' Be sure to open this in a scripting text editor.
1. Trim and align raw reads, remove PCR duplicates >  get library compositions, RPF length distributions, bedgraphs for viewing on IGV, plot RPF length distributions
2. Get periodicity for each length of read > 5' and 3' frame counts files. This step is optional and is not used downstream, but it is useful for deciding read lengths and offsets to include and use in A site transformation
3. Make bed files for Vplots  >  5' and 3' files for input to ScatterPlotsLengths.R to visualize read lengths along genes. This can be used along with step 2 to determine A site offset.
4. A-site transformation  >  get periodicity and coverage, A-site bedgraphs for viewing on IGV and downstream analyses
5. Normalize read counts to reads per million mapped (rpm) > Standalone script to rpm-normalize, though this function is also included in some steps downstream. 
6. Combine data from replicates > sum all data, and also rpm-normalize the combined file
7. Get enrichment over control sample > Use rpm-normalized files (summed or individual) to calculate enrichment in a sliding window. Window size and coverage threshold are taken as inputs. Output includes both log2 and non-log2 transformed values.
8. Various plotting scripts in R to make visualizations and quantify read coverage



# Data availability and manuscript

Fastq files are deposited in the GEO database under the accession number GSEXXXXXX. The link to our full manuscript will be provided here upon publication.




To install:  
Download zipped repository and unzip into directory with *fastq.gz files.
Install R packages as needed, listed in each R script.  
To run:
Open yMitoRP_pipelineCommands_README.txt in a scripting text editor and follow examples, adjusting each script (paths, etc) as necessary.



 

