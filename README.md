# IsoQMap

![PyPI](https://img.shields.io/pypi/v/isoqmap)
![Build](https://img.shields.io/badge/build-passing-brightgreen)


**IsoQMap** is an automated pipeline for isoform expression quantification from RNA-seq data and subsequent isoform-level QTL (isoQTL) mapping. It integrates two powerful tools:

- **[XAEM](https://github.com/WenjiangDeng/XAEM)** – a robust method for isoform expression estimation across samples ([paper](https://academic.oup.com/bioinformatics/article/36/3/805/5545974), [website](https://www.meb.ki.se/sites/biostatwiki/xaem)).
- **[OSCA](https://yanglab.westlake.edu.cn/software/osca/)** – for genetic mapping of isoforms and genes using multi-omics data ([paper:OSCA](https://pubmed.ncbi.nlm.nih.gov/31138268/), [paper:THISTLE](https://www.nature.com/articles/s41588-022-01154-4), [website](https://yanglab.westlake.edu.cn/software/osca/))

---

## 📦 Prerequisites

- Python ≥ 3.8
- R ≥ 3.6.1

---

## 🛠️ Installation

### Using `conda` prepare prerequisites (Recommended)
```bash
conda create -n IsoQMap python=3.8 r-base=4.1.2 r-essentials
conda activate IsoQMap
conda install -c conda-forge r-foreach r-doparallel
```

### Quick installation
```bash
pip install isoqmap
# (For China Mainland)
pip install isoqmap -i https://pypi.tuna.tsinghua.edu.cn/simple
```

---

## 🚀 Quick Start

```bash
isoqmap --help
```

---

## 📁 Example

A working example command is provided in the `examples/` directory:

```bash
cd gitpath/to/isoqmap/examples
sh Demo1.run.isoqmap_pipeline.sh
```

---

## 🔬 Isoform Expression Quantification (`isoqmap isoquan`)

### 🔹 Prepare Input file

Prepare a tab-delimited file for FASTQ file (e.g., `infastq_lst.tsv`) with four columns:

```
SampleName   SourceName   FASTQ_R1   FASTQ_R2
```

#### Example: Single Batch
```
sample4   S0007   S0007_1.fq.gz   S0007_2.fq.gz
sample5   S0008   S0008_1.fq.gz   S0008_2.fq.gz
```

#### Example: Multiple Batches

```
sample1   S0001   S0001_1.fq.gz   S0001_2.fq.gz
sample1   S0002   S0002_1.fq.gz   S0002_2.fq.gz
sample2   S0003   S0003_1.fq.gz   S0003_2.fq.gz
sample2   S0004   S0004_1.fq.gz   S0004_2.fq.gz
```

### 🔹 Run isoquan 

```bash
isoqmap isoquan -i /path/to/infastq_lst.tsv
```

#### Optional:

- Specify a reference:
  ```bash
  --ref gencode_38
  ```
- Provide a custom config:
  ```bash
  -c /path/to/config.ini
  ```
- Specify a out directory
  ```bash
   -o path/to/workdir
  ```
---

## 🧬 Isoform and Gene QTL Mapping (`isoqmap isoqtl`)

### Pipline quick start (precess -> qtl calling -> format)
```bash
isoqmap isoqtl pipeline -i   /path/to/XAEM_isoform_expression_tpm.tsv.gz  --bfile  /path/to/genotype_for_isoqmap --covariates QTL_covariate.tsv 
```
#### bfile demo
```text
genotype_for_isoqmap.bed  genotype_for_isoqmap.bim  genotype_for_isoqmap.fam
```
#### covariates demo
| ID | S367402 | S367403 | ... | S367404 |
|----|---------|---------|-----|---------|
| Sex | 2       | 1       | ... | 2       |
| age | 64      | 54      | ... | 65      |
| ... | ...     | ...     | ... | ...     |
| PC1 | -0.1037 | 0.0118  | ... | 0.0112  |
| PC2 | -0.0173 | -0.003  | ... | -0.008  |
| PC3 | 0.0017  | 0.0008  | ... | 0.0308  |


#### Optional:
- Specify a reference:
  ```bash
  --ref gencode_38
  ```
- Provide a custom config:
  ```bash
  -c /path/to/config.ini
  ```
- Specify a out directory
  ```bash
   -o path/to/workdir
  ```
- Specify a prefix
  ```bash
   -outprefix qtl_job
  ``` 
---

### Run by each model
#### Step 1: Preprocess input files for QTL mapping
```bash
outdir=workdir
isoqmap isoqtl preprocess -i path/to/XAEM_isoform_expression_tpm.tsv.gz --isoform-ratio --ref gencode_38 --covariates QTL_covariate.tsv --outdir $outdir 
```
This step involves transforming isoform expression data into isoform ratios, applying normalization, adjusting for covariates, and generating the input BOD file for downstream QTL mapping.

#### Step 2: Run QTL mapping (eQTL / isoQTL / irQTL)

##### eQTL
```bash
isoqmap isoqtl call --bfile genotype/test_for_isoqmap --befile $outdir/BOD_files/IsoQ.gene_abundance --mode eqtl --outdir $outdir/QTL_results  --run   
```                                               
##### isoQTL
```bash
isoqmap isoqtl call --bfile genotype/test_for_isoqmap --befile $outdir/BOD_files/IsoQ.isoform_abundance --mode sqtl --outdir $outdir/QTL_results  --run     
```                                               
##### irQTL
```bash
isoqmap isoqtl call --bfile genotype/test_for_isoqmap --befile $outdir/BOD_files/IsoQ.isoform_splice_ratio --mode sqtl --outdir $outdir/QTL_results  --run
```

#### Step 3: Format QTL results
#### isoQTL and irQTL
```bash
isoqmap isoqtl format --infile "$outdir/QTL_results/osca_qtl_job.*.sqtl_10_*_isoform_eQTL_effect.txt" --mode sqtl --ref gencode_38      
```
##### 
```bash
isoqmap isoqtl format --infile "$outdir/QTL_results/osca_qtl_job*eqtl_10_*.besd" --mode eqtl --ref gencode_38  
```
This step formats the results for downstream Mendelian Randomization (MR), Colocalization (coloc), or other integrative analyses.

---

## 📬 Feedback

For issues, bug reports, or feature requests, please open an issue or submit a pull request.

---

## 📄 License

MIT License