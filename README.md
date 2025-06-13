# IsoQMap

**IsoQMap** is an automated pipeline for isoform expression quantification from RNA-seq data and subsequent isoform-level QTL (isoQTL) mapping. It integrates two powerful tools:

- **[XAEM](https://github.com/WenjiangDeng/XAEM)** – a robust method for isoform expression estimation across samples ([paper](https://academic.oup.com/bioinformatics/article/36/3/805/5545974), [website](https://www.meb.ki.se/sites/biostatwiki/xaem)).
- **[OSCA](https://yanglab.westlake.edu.cn/software/osca/)** – for genetic mapping of isoforms and genes using multi-omics data.

---

## 📦 Prerequisites

- Python ≥ 3.8
- R ≥ 3.6.1

---

## 🛠️ Installation

### Using `conda`:

```bash
conda create -n IsoQMap python=3.8 r-base=4.1.2 r-essentials
conda activate IsoQMap
conda install -c conda-forge r-foreach r-doparallel
pip install isoqmap
```

---

## 🚀 Quick Start

```bash
isoqmap --help
```

---

## 📁 Example

A working example is provided in the `Example/` directory:

```bash
cd /path/to/isoqmap/Example
sh run_example.sh
```

---

## 🔬 Isoform Expression Quantification (`isoqmap isoquan`)

### 🔹 Input Format

Prepare a tab-delimited file (e.g., `infastq_lst.tsv`) with four columns:

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

### 🔹 Run XAEM

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

### Step 1: Convert Expression BOD file
```bash
isoqmap isoqtl preprocess
```

---

### Step 2: Run QTL Mapping

```bash
isoqmap isoqtl call 
```
### Step 3: Format for QTL results for MR and Coloc
```bash
isoqmap isoqtl format 
```
---

## 📬 Feedback

For issues, bug reports, or feature requests, please open an issue or submit a pull request.

---

## 📄 License

MIT License