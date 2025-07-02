outdir=isoqmap_analysis
#isoform quanlification
isoqmap isoquan -i fq.list \
                -c config.ini \
                -o $outdir \
                --ref gencode_38

#isoform qtl mapping
## data propcoess
isoqmap isoqtl preprocess -i $outdir/results/XAEM_isoform_expression_tpm.tsv.gz \
                          --isoform-ratio --ref gencode_38 \
                          --covariates QTL_covariate.tsv \
                          --outdir $outdir
## qtl calling
### eQTL
isoqmap isoqtl call --bfile genotype/test_for_isoqmap \
                    --befile $outdir/BOD_files/IsoQ.gene_abundance \
                    --mode eqtl \
                    --outdir $outdir/QTL_results  \
                    --run
#### isoQTL
isoqmap isoqtl call --bfile genotype/test_for_isoqmap \
                    --befile $outdir/BOD_files/IsoQ.isoform_abundance \
                    --mode sqtl \
                    --outdir $outdir/QTL_results  \
                    --run
#### irQTL
isoqmap isoqtl call --bfile genotype/test_for_isoqmap \
                    --befile $outdir/BOD_files/IsoQ.isoform_splice_ratio \
                    --mode sqtl --outdir $outdir/QTL_results  \
                    --run
### qtl format
#### isoQTL and irQTL
isoqmap isoqtl format --infile "$outdir/QTL_results/osca_qtl_job.*.sqtl_10_*_isoform_eQTL_effect.txt" \
                      --mode sqtl \
                      --ref gencode_38
#### eQTL
isoqmap isoqtl format --infile "$outdir/QTL_results/osca_qtl_job*eqtl_10_*.besd" \
                      --mode eqtl \
                      --ref gencode_38

