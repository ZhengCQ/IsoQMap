# SpReD-GIT 
Splicing-regulatory Driver Genes Identification Tool

This project describes a tool to identify the splicing-regulatory driver genes based on isoform-level expressions quanlification such as XAEM, RSEM

## Prerequisites
```
Python (recommended version >= 3.7)
```

## Installation

#### Step 1: Setup tools
```
git clone https://github.com/ZhengCQ/SpReD.git
cd SpReD
pip install -e ./
```


## Demo
```
spred run-all -m matrix.count.input.tsv -e group.tsv -g Group_name -c1 Case_name -c2 Control_name -o outdir_name --species mouse/human
```# IsoQMap
