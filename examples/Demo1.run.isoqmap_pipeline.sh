outdir=isoqmap_analysis
#isoform quanlification
isoqmap isoquan -i fq.list \
                -c config.ini \
                -o $outdir \
                --ref gencode_38

#isoform qtl mapping
isoqmap isoqtl pipeline -i  $outdir/results/XAEM_isoform_expression_tpm.tsv.gz \
                        --bfile genotype/test_for_isoqmap \
                        --ref gencode_38 \
                        --covariates QTL_covariate.tsv \
                        --outdir $outdir
                        