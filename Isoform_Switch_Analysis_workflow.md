## Isoform Switch Analysis Workflow

### Project's goal: detect isoform switching events between conditions

### Workflow steps:

1. Initial quality control (FastQC https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. Adapter and quality trimming, if deemed necessary (Trim Galore! https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
3. Isoform quantification (RSEM https://deweylab.github.io/RSEM/)
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

### 3. Isoform quantification with RSEM.

RSEM paper: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323

RSEM tutorial: https://github.com/bli25broad/RSEM_tutorial

1) First we will prepare the reference. RSEM works with a set of transcripts, that can be prepared using the reference genome and a GTF file with exon coordinates.

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

To prepare the transcript reference run the following command:

```
mkdir ref
rsem-prepare-reference --gtf <path/to/gtf> --bowtie2 </path/to/genome.fa> ref/

```
This command assumes that we will use bowtie2 as internal aligner. The results will be saved in ref/ directory.

Once the refernce is prepared, we can run RSEM:

```
rsem-calculate-expression -p 24 --bowtie2 --estimate-rspd --append-names --output-genome-bam <FASTQ> ../ref/human_0 <SAMPLE_NAME>

```

Iterate over all of the fastq files in the directory and run RSEM isoform quantification.

```
#! /bin/bash

dir=$1

for file in $dir/*.fastq
do
    filename=`basename $file`
    samplename=${filename%.*}
    rsem-calculate-expression -p 24 --bowtie2 --estimate-rspd --append-names --output-genome-bam $file ../ref/human_0  $samplename

done

```
NOTE! In one of the conda environments I encountered the a bowtie2 installation problem breaking RSEM run.
The solution to the problem is described in the following Biostars post:
https://www.biostars.org/p/494922/

See the post content below:

```

$ conda create -n bttest -c bioconda bowtie2
$ conda activate bttest
$ bowtie2
~/.conda/envs/bttest/bin/bowtie2-align-s: error while loading shared libraries: libtbb.so.2: cannot open shared object file: No such file or directory
(ERR): Description of arguments failed!
Exiting now ...

```

The solution is to downgrade the tbb version

```
$ conda install tbb=2020.2
$ bowtie2 -h

```




















