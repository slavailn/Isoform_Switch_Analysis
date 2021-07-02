## Isoform Switch Analysis Workflow

### Software installation
Use conda to install the following packages:

```
conda install -c bioconda fastqc
conda install -c bioconda trim_galore
conda install -c bioconda gffread
conda install -c bioconda kallisto

```

### Project's goal: detect isoform switching events between conditions

### Workflow steps:

1. Initial quality control (FastQC https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. Adapter and quality trimming, if deemed necessary (Trim Galore! https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
3. Isoform quantification (Kallisto https://pachterlab.github.io/kallisto/)
4. Detection of alternative splicing and isoform switch events (IsoformSwicthAnalyzeR https://www.bioconductor.org/packages/release/bioc/html/IsoformSwitchAnalyzeR.html)

### 1. Initial quality control

Run quality analysis using FastQC software.

Assuming that **fastqc** is in the system's path, change to the directory containing FASTQ files and run:

```

fastqc *.fastqc

```

To run FastQC analysis in parallel:

```

parallel fastqc {} ::: *.fastq

``` 

Examine the results by looking into html report files. They can be opened in any browser.
The most relevant sections of the quality report are basic statistics telling the total number of the reads; distributions of base qualities; per base sequence content;
GC content; over-represented sequences. The presence of adapters used in the library preparation process is evident from over-represented sequences report section. We will skip adapter trimming if no adapters are present.

FastQC tutorial: https://www.youtube.com/watch?v=bz93ReOv87Y
Site with many examples of sequencing QC problems: https://sequencing.qcfail.com/

### 2. Adapter trimming (if necessary)

```
trim_galore --illumina <FASTQ> # for Illumina Truseq adapters

```

See Trim Galore help for other options.

### 3. Isoform quantification with Kallisto.

Kallisto paper: https://www.nature.com/articles/nbt.3519

Kallisto tutorial: https://pachterlab.github.io/kallisto/starting

1) First we will prepare the reference. Kallisto works with a set of transcripts, that can be prepared using the reference genome and a GTF file with exon coordinates.

An easy way to obtaine an annotated genome is to download one from Illumina iGENOME repository - https://support.illumina.com/sequencing/sequencing_software/igenome.html

For example, get an Ensembl assembly of a human genome with wget

```
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/Ensembl/GRCh37/Homo_sapiens_Ensembl_GRCh37.tar.gz

```
Untar and unzip the genome tarball.

```

tar -xzvf Homo_sapiens_Ensembl_GRCh37.tar.gz

```

The GTF file is located in Annotation/Genes/ directory and the whole genome fasta file is found in the Sequence/ sub-folder.

To prepare the transcript reference we will first extract the transcript sequences using gffread utility: http://ccb.jhu.edu/software/stringtie/gff.shtml

```
mkdir ref
gffread genes.gtf -g <Genome.fasta> -w <transcripts.fasta>
# Now create kallisto index, we will use default 31 nt k-mer length.
kallisto index -i <transcripts.idx> <transcripts.fa>

```

Once the reference index is prepared, we can run quantify transcript abundances with kalliso:

```
kallisto quant -i </path/to/index> \ # path to index file (.idx) 
               -o <out_dir> \ path to output directory 
			   -b 100 \ # number of bootstrap samples
			   -bias \ # perform sequence bias correction
 			   --single \ # quantify single end reads
			   -l 300 \ # estimated average fragment length
			   -s 20 \ # estimated standard deviation fragment length
			   <FASTQ> # sample fastq

```

Iterate over all of the fastq files in the directory and run Kallisto isoform quantification.

```
#! /bin/bash

# Provide path to directory with fastq files as the first argument
# and path to kallisto index (.idx) file as the second argument.

dir=$1 # path to directory with fastq files 
idx=$2 # path kallisto index file

for file in $dir/*.fastq
do
    filename=`basename $file`
    samplename=${filename%.*}
    echo -e "#------ Processing $samplename ------- #"
    kallisto quant -t 20 -i $2 -o $samplename -b 100 -bias --single -l 300 -s 20 $file

done

echo -e "ALL FILES FINISHED!"

```

Kallisto will generate a folder with the output files for each of the samples. 
The folder contains 2 files with transcript abundances and a run info json file.

The tab delimited .tsv file has colums that correspond to transcrip ids, transcript length (actual and estimated), estimated counts and transcripts per million (TPM).

```
target_id       length  eff_length      est_counts      tpm
ENST00000456328 1657    1358    0       0
ENST00000515242 1653    1354    0       0
ENST00000518655 1483    1184    0       0
ENST00000450305 632     333     0       0
ENST00000541675 1416    1117    383.017 29.2701
ENST00000423562 1669    1370    96.199  5.99389
ENST00000438504 1783    1484    0       0
ENST00000488147 1351    1052    0       0
ENST00000538476 1583    1284    0       0
```

### 4. Detect, visualize, and annotate isoform switching events with IsoSwitchAnalyzeR

Follow the R workflow as shown below. 
The transcript abundances that serve as input into IsoformSwicthAnalyzeR were obtained using kallisto and saved in **kallisto_results** directory.

```
library(IsoformSwitchAnalyzeR)

setwd("<path/to/project_dir>")

### This example is based on comparing CT_18 and CT_DMSO experimental groups.

### Import Kallisto isoform data
setwd("kallisto_results/")
list.files()
kallisto_files <- list.files(".", pattern = "tsv", recursive = T)
kallisto_files

# [1] "CT_18_1/abundance.tsv"   "CT_18_2/abundance.tsv"   "CT_18_3/abundance.tsv"  
# [4] "CT_DMSO_1/abundance.tsv" "CT_DMSO_2/abundance.tsv" "CT_DMSO_3/abundance.tsv"

kallistoQuant <- importIsoformExpression(
  sampleVector = kallisto_files,
  addIsofomIdAsColumn = TRUE,
  interLibNormTxPM = TRUE,
  normalizationMethod = 'TMM',
  showProgress = TRUE
  
)

kallistoQuant
setwd("../")
list.files()
save(kallistoQuant, file="RData_objects/kallistoQuant.RData")

### Create design matrix
colnames(kallistoQuant$abundance)
conditions <- c("CT_18", "CT_18", "CT_18", "CT_DMSO", 
                "CT_DMSO", "CT_DMSO") 
conditions

myDesign <- data.frame(
  sampleID = colnames(kallistoQuant$abundance)[-1],
  condition = conditions)
myDesign

### Create isoformSwitchList object
aSwitchList <- importRdata(
  isoformCountMatrix   = kallistoQuant$counts, # counts for statistical analysis
  isoformRepExpression = kallistoQuant$abundance, # abundances for effect size calculations 
  designMatrix         = myDesign, # design matrix
  isoformExonAnnoation = "transcript_annot/genes.gtf", # gtf file used to generate transcripts
  isoformNtFasta       = "transcript_annot/transcripts.fa", # fasta file with transcript sequences
  showProgress = T
)

aSwitchList
save(kallistoQuant, file="RData_objects/kallistoQuant.RData")

# run the first part of the isoform switch analysis workflow which filters for 
# non-expressed genes/isoforms, identifies isoform switches, 
# annotates open reading frames (ORF), switches with and extracts both 
# the nucleotide and peptide (amino acid) sequences and output them as 
# two separate fasta files 

## Prefilter genes with low expression and remove single isoform genes
aSwitchListFiltered <- preFilter(
  switchAnalyzeRlist = aSwitchList,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE
)

IsoSwitchList <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist   = aSwitchListFiltered,
  outputSequences = TRUE, # output sequences for annotation 
  prepareForWebServers = TRUE  # for sequence analysis on web servers
)

# Save the results summary
switchSummary <- extractSwitchSummary( IsoSwitchList )
switchSummary
write.table(switchSummary, file = "switch_summary.csv",
            sep = ",", col.names = T, row.names = F)
save(IsoSwitchList, file="IsoSwitchList_part1_compleleted.RData")

```

























