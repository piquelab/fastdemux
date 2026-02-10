# fastdemux

Preprint available at 

# Installation

```
git clone https://github.com/piquelab/fastdemux.git
cd fastdemux
git clone https://github.com/attractivechaos/klib.git
mkdir build
cd build
cmake ../
make
make install
```
Note that we require [htslib](https://github.com/samtools/htslib) to be installed. You may need to adjust the library location in the [CMakeLists.txt](CMakeLists.txt). It also uses OpenMP and Zlib that should be standard on most systems. You also need [klib](https://github.com/samtools/klib) which can be cloned within the repository. We also provide a Google Colab example showing how to install and use [testColab.ipynb](testColab.ipynb). 


# Running fastdemux

`fastdemux` uses a DLDA approach to assign single cells to donors using genotype information. This section describes how to run `fastdemux` on Cell Ranger output using a donor genotype VCF.

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
  For example, if the BAM uses `chr1`…`chr22`, the VCF should also use `chr1`…`chr22` with the same chromosome ordering. Mismatches in chromosome annotation or ordering can lead to failed runs or incorrect results. It is important to independently check that the vcf and bam files use the same reference genome.
- **Recommended VCF preprocessing:** filter to **biallelic SNPs** and prefer sites with **≥10× read coverage** to improve robustness of demultiplexing.

Note that it is preferrable that a large amount of SNPs with coverage be imputed. Dosage DS values are preferred but GP or GT values can be used if not available. You can use GENCOVE, GLIMPSE2, or TopMed Imputation server, to generage an imputed vcf file from low coverage sequencing data or from genotyping microarrays. We have repos ([scRNAseq repo](https://github.com/piquelab/counts_cellranger), [scATACseq repo](https://github.com/piquelab/counts_cellranger_atac)) with templates on how we use cellranger and fastdemux in a production environment at the Wayne State University High Performance Computing Grid. The scripts show how to also filter the genotype VCF file and check that bam and vcf are consistently ordered. 


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

## Benchmarking

We benchmarked `fastdemux` against commonly used genotype-based demultiplexing tools (demuxlet, vireo, and demuxalot) across varying sequencing depths and genotype SNP coverage thresholds. Performance was evaluated in terms of donor assignment error rate, runtime, and peak memory usage.

Across all tools, donor assignment error rates decreased as sequencing depth increased. At very low read fractions (1–5%), all methods performed similarly. At moderate to high read fractions (≥30%), `fastdemux` consistently achieved the lowest total-droplet error rates, indicating improved robustness as sequencing depth increases. In addition to improved accuracy, `fastdemux` showed substantially lower runtime and memory usage than demuxlet and vireo, with memory usage remaining near constant across read depths.

![Benchmarking across read depth](https://raw.githubusercontent.com/piquelab/fastdemux_bench/main/figures/fig2new.png)
![Benchmarking across read depth](fig2new.png)

When varying the minimum SNP coverage threshold used to filter the genotype VCF, all tools exhibited increasing error rates as fewer SNPs were retained. Across all thresholds, `fastdemux` consistently achieved lower error rates than demuxlet, vireo, and demuxalot. Performance differences were most pronounced at lower SNP coverage thresholds (G49 and G9), where `fastdemux` retained the lowest error rate while maintaining fast runtimes and minimal memory usage. In contrast, demuxlet and vireo showed sharp increases in runtime and memory as SNP density increased.

![Benchmarking across SNP coverage thresholds](https://raw.githubusercontent.com/piquelab/fastdemux_bench/main/figures/fig3new.png)
![Benchmarking across SNP coverage thresholds](fig3new.png)

Overall, these benchmarks demonstrate that `fastdemux` effectively leverages large numbers of lower-coverage SNPs to achieve accurate donor assignment while remaining computationally efficient, making it well suited for large-scale and low-coverage single-cell datasets, including scATAC-seq.

