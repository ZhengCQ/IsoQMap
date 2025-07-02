outdir=isoqmap_analysis_pipeline
#isoform quanlification
isoqmap isoquan -i fq.list \
                -c config.ini \
                -o $outdir \
                --ref gencode_38

#isoform qtl mappint
isoqmap isoqtl pipeline -i  isoqmap_analysis/results/XAEM_isoform_expression_tpm.tsv.gz \
                        --bfile genotype/test_for_isoqmap \
                        --ref gencode_38 \
                        --covariates QTL_covariate.tsv \
                        --outdir $outdir
