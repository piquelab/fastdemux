# fastdemux

# Installation

# Running fastdemux

`fastdemux` uses a DLDA approach to assign single cells to donors using genotype information. This section describes how to run `fastdemux` on Cell Ranger outputs using a donor genotype VCF.

## Required inputs

For each **library**, you need:

- **Aligned BAM file**:  
  `possorted_genome_bam.bam`
- **Barcode list**:  
  `raw_feature_bc_matrix/barcodes.tsv.gz`  
  (or `filtered_feature_bc_matrix/barcodes.tsv.gz`)
- **Genotype VCF** containing donor genotypes:  
  `*.vcf.gz` (bgzipped and index file)

### Important requirements and recommendations

- **Chromosome naming and ordering must match between the BAM and VCF.**  
  For example, if the BAM uses `chr1`…`chr22`, the VCF should also use `chr1`…`chr22` with the same chromosome ordering. Mismatches in chromosome annotation or ordering can lead to failed runs or incorrect results.
- **Recommended VCF preprocessing:** filter to **biallelic SNPs** and prefer sites with **≥10× read coverage** to improve robustness of demultiplexing.

## Example directory structure

```text
project/
├── genotypes/
│   ├── combined.vcf.gz
│   └── combined.vcf.gz.tbi
├── libraryA/
│   ├── possorted_genome_bam.bam
│   └── raw_feature_bc_matrix/
│       └── barcodes.tsv.gz
├── libraryB/
│   ├── possorted_genome_bam.bam
│   └── raw_feature_bc_matrix/
│       └── barcodes.tsv.gz
└── fdout/
```

## Basic usage

`fastdemux` takes four positional arguments:

1. BAM file  
2. Genotype VCF  
3. Barcode file  
4. Output prefix  

Example command:

```bash
fastdemux -t 2 \
  libraryA/possorted_genome_bam.bam \
  genotypes/combined.vcf.gz \
  libraryA/raw_feature_bc_matrix/barcodes.tsv.gz \
  fdout/libraryA.fdout.raw
```

### Parameters

- `-t 2`  
  Number of threads to use.
- **Output prefix**  
  `fastdemux` will generate output files in the `fdout` directory using this specified prefix.

## Output files and interpretation

`fastdemux` produces four output files per library that summarize donor similarity, donor assignment, assignment confidence, and model scores.

### 1. Correlation file (`*.fdout.raw.corr.txt.gz`)
This file reports pairwise correlations between inferred donor genotype profiles and is primarily used as a diagnostic to assess donor similarity and genotype distinguishability.

**Columns:**
- `Sample1`: donor identifier for the first sample in the pair
- `Sample2`: donor identifier for the second sample in the pair
- `Correlation`: correlation coefficient between the genotype signals of the two donors

### 2. Raw assignment file (`*.fdout.raw.info.txt.gz`)
This file contains donor assignment information for all barcodes in the library, including low-information droplets, and is intended for quality control and diagnostic inspection.

### 3. Filtered assignment file (`*.fdout.raw.info.2nd.txt.gz`)
This file contains a subset of barcodes with sufficient allelic information for confident demultiplexing and is used for downstream analyses.

**Columns:**
- `BARCODE`: cell barcode identifier
- `bcnum`: internal numeric index assigned to the barcode
- `Nsnp`: number of informative SNPs observed for the barcode
- `Numi`: number of UMIs for the assigned cell
- `DropType`: droplet classification; 1=singlet, 2=doublet, 3=triplet, etc. 
- `KletScore`: likelihood-based score summarizing genotype consistency across donors
- `BestScore`: assignment score for the most likely donor
- `BestSample`: donor assigned to the barcode
- `SecondBestScore`: assignment score for the second most likely donor
- `SecondBestSample`: donor with the second-highest assignment score
- `ThirdBestScore`: assignment score for the third most likely donor
- `ThirdBestSample`: donor with the third highest assignment score

### 4. DLDA score file (`*.fdout.raw.dlda.txt.gz`)
This file contains donor-specific discriminant scores for each barcode and represents the core model output used internally to determine donor assignments.

**Columns:**
- `BARCODE`: cell barcode identifier
- one column per donor (e.g. `NA18486`, `NA18489`, …): discriminant score for assigning the barcode to the corresponding donor, with higher values indicating stronger support

### Recommended usage

- Use `*.fdout.raw.info.2nd.txt.gz` for downstream analyses requiring high-confidence donor assignments.
- Use `*.fdout.raw.info.txt.gz` for quality control, diagnostics, and inspection of low-information barcodes.
- Use `*.corr.txt.gz` and `*.dlda.txt.gz` primarily for diagnostic, benchmarking, and methodological evaluation.
