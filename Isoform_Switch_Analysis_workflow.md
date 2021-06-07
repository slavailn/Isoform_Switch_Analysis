## Isoform Switch Aanalysis Workflow

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

### 2. Adapter trimming (if necessary)

```
trim_galore --illumina <FASTQ> # for Illumina Truseq adapters

```

See Trim Galore help for other options.

### 3. Isoform quantification with RSEM.










