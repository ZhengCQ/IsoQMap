# IsoQMap
This project describes an automated pipeline to quantify isoform expressions for the RNA sequencing data and then for isoform QTL mapping. The protocol uses [XAEM](https://github.com/WenjiangDeng/XAEM), a powerful method for isoform expression estimation across multiple samples. You also can find the detailed information about XAEM [website](https://www.meb.ki.se/sites/biostatwiki/xaem) and in the published paper in [Bioinformatics](https://academic.oup.com/bioinformatics/article/36/3/805/5545974).

## Prerequisites
```
R (recommended version >= 3.5.1)
Python (recommended version >= 3.7)
```

## Development


## Installation
### conda
```
conda create -n IsoQMap python=3.8  r-base=4.1.2 r-essentials
conda activate IsoQMap
conda install -c conda-forge r-foreach r-doparallel
pip install isoqmap
``` 

## qucik start
```
isoqmap --help
```

## Example
An example is prepared in the project **Example** folder, executable as
```
cd /path/to/iGTEx_XAEM/Example
sh run_example.sh 
```

## isoqmap
#### Input files
In `/path/to/workdir`, create a file `/path/to/workdir/infastq_lst.tsv` listing the FASTQ input files. The file is a tab-delimited text file with 4 columns: `Sample name`, `Source name`, `FASTQ file name for paired-end read 1`, and `FASTQ file name for paired-end read 2`. `Source name` indicates the batch or sequencing library of the sample, so that the same sample may correspond to more than one sources. A standard example, where each sample has only a single batch or multiple batches, is given as `/path/to/isoqmap/Example/infastq_lst.tsv`:

##### single batch
```
sample4 S0007   S0007_1.fg.gz   S0007_2.fg.gz
sample5 S0008   S0008_1.fg.gz   S0008_2.fg.gz
```
##### multiple batches
```
sample1 S0001   S0001_1.fg.gz   S0001_2.fg.gz
sample1 S0002   S0002_1.fg.gz   S0002_2.fg.gz
sample2 S0003   S0003_1.fg.gz   S0003_2.fg.gz
sample2 S0004   S0004_1.fg.gz   S0004_2.fg.gz
sample3 S0005   S0005_1.fg.gz   S0005_2.fg.gz
sample3 S0006   S0006_1.fg.gz   S0006_2.fg.gz
```


#### Run isoform quanlification 
XAEM can be easily run with:
```
isoqmap isoquan -i /path/to/Tissue1/infastq_lst.tsv
```
(Optional) To specify a particular reference with gencode hg38, use:
```
--ref gencode_38
```
(Optional) To specify a particular output directory, use:
```
-o /path/to/Tissue1_output_directory
```
(Optional) Further customized configuration of XAEM can be setup by:
```
-c /path/to/Tissue1_config.ini
```
An example of the `config.ini` file can be found in `/path/to/isoqmap/config.ini/`.


## isoqtl

### preprocess
Prepare for QTL mapping fot eQTL/isoQTL/irQTL, from tpm expression matrix to BOD file
```
isoqmap isoqtl preprocess -i path/to/XAEM_isoform_expression_tpm.tsv.gz --isoform-ratio --ref gencode_38 --covariates /path/to/QTL_covariates.tsv
```
### call

### format
#### sqtl
```
isoqmap isoqtl format --infile 'workdir/QTL_results/osca_qtl.*.sqtl_10_*_isoform_eQTL_effect.txt.gz' --mode sqtl --ref gencode_38 --id2rs-file path/to/anno_into.tsv.gz
```
#### eqtl
isoqmap isoqtl format --infile 'workdir/QTL_results/osca_qtl.gene_abundance.eqtl_10_*.besd' --mode eqtl --ref gencode_38 --id2rs-file path/to/anno_into.tsv.gz