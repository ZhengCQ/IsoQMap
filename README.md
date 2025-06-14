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

A working example command is provided in the `Example/` directory:

```bash
cd /path/to/isoqmap/Example
sh run_example.sh
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

---

## 🧬 Isoform and Gene QTL Mapping (`isoqmap isoqtl`)

### Step 1: Preprocess input files for QTL mapping
```bash
isoqmap isoqtl preprocess
```
This includes converting isoform expression to isoform ratio, normalization, preadjust covaries, and preparing input (BOD file) for QTL mapping

### Step 2: Run QTL mapping (eQTL / isoQTL / irQTL)
```bash
isoqmap isoqtl call 
```
You can specify QTL types and models via CLI options.


### Step 3: Format QTL results
```bash
isoqmap isoqtl format 
```
This step formats the results for downstream Mendelian Randomization (MR), Colocalization (coloc), or other integrative analyses.

---

## 📬 Feedback

For issues, bug reports, or feature requests, please open an issue or submit a pull request.

---

## 📄 License

MIT License
