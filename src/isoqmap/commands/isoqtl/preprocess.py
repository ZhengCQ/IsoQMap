import os
import sys
import click
import pandas as pd
import numpy as np
import datetime
import scipy.stats as stats
import statsmodels.api as sm
from ...tools import pathfinder,common
import logging
import subprocess

from ...tools.downloader import download_reference, download_osca



# Configure logging
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)
binfinder = pathfinder.BinPathFinder('isomap')


def check_input_files(isoform_file, cov_file, refdb):
    if isoform_file.endswith('RData'):
        script_path = binfinder.find('tools/isoform_rdata2exp.R')
        os.system(f'Rscript {script_path} inRdata={isoform_file}')
        isoform_file = isoform_file.replace(".RData", "_tpm.tsv")
        logger.info(f'Converted RData to: {isoform_file}')
    
    ## covariate
    logger.info(f'Reading covariate file {cov_file} and checking covariate samples')
    df_cov = pd.read_csv(cov_file, sep='\t', index_col=0)
    df_cov.columns = [i.replace(' ', '') for i in df_cov.columns]
    logger.info(f'There are {df_cov.shape[1]} Samples and {df_cov.shape[0]} items for adjust in covariate file.')
    ### Read expression file for columns
    logger.info(f'Prereading expression file {isoform_file}...')
    df_exp_head = pd.read_csv(isoform_file, sep='\t', index_col=0,nrows=1)
    logger.info(f'There are {df_exp_head.shape[1]} Samples in expression file.')

    
    common_samples = df_exp_head.columns.intersection(df_cov.columns)
    match_ratio = len(common_samples) / df_cov.shape[1]

    if match_ratio < 0.8:
        logger.error(f'Only {match_ratio:.2%} of covariate samples are found in isoform expression matrix. Aborting.')
        sys.exit(1)
    elif match_ratio < 1.0:
        logger.warning(f'{len(common_samples)} out of {df_exp_head.shape[1]} samples matched. Filtering unmatched samples.')
    elif match_ratio == 1:
        logger.warning(f'{len(common_samples)} out of {df_exp_head.shape[1]} samples matched. Keep going.')

    
    df_cov = df_cov[common_samples]
    ## expression file
    logger.info(f'Re reading expression file {isoform_file}...')
    df_exp = pd.read_csv(isoform_file, sep='\t', index_col=0)[common_samples]
    logger.info(f'There are {df_exp.shape[1]} Samples and {df_exp.shape[0]} isoform transcripts in expression file.')

    
    ## annotation file
    gene_info_fi = str(binfinder.find(f'./resources/ref/{refdb}/transcript_gene_info.tsv.gz'))    
        
    if not common.check_file_exists(
        gene_info_fi,
        file_description=f"Gene annotaion file {gene_info_fi}",
        logger=logger,
        exit_on_error=False
    ):
        print(f"Gene annotaion file not found or unreadable. Trying to download for {refdb}...")
        download_reference(refdb, ['geneinfo'])
    
        gene_info_fi = str(binfinder.find(f'./resources/ref/{refdb}/transcript_gene_info.tsv.gz'))    
 

    logger.info(f'Reading annotation file {gene_info_fi}...')
    try:
        df_anno = pd.read_csv(gene_info_fi, sep='\t',index_col=0)
    except:
        logging.ERROR(f"Error for eading annotation file {gene_info_fi}")

    common_transcript = df_exp.index.intersection(df_anno.index)
    logger.warning(f'{len(common_transcript)} out of {df_exp.shape[0]} transcript matched. Filtering unmatched trascripts.')
    
    df_exp = df_exp.loc[common_transcript]

    return df_exp, df_anno, df_cov


def filtered_isoform(df_exp, df_anno, tpm_threshold=0.1, sample_threshold_ratio=0.2):
    threshold_count = df_exp.shape[1] * sample_threshold_ratio
    low_exp_mask = (df_exp >= tpm_threshold).sum(axis=1) <= threshold_count
    logger.info(f'{low_exp_mask.sum()} isoforms excluded (TPM < {tpm_threshold} in > {100*(1-sample_threshold_ratio):.0f}% samples)')

    df_exp = df_exp[~low_exp_mask]

    df_exp = df_anno[['gene_id']].merge(df_exp, left_index=True, right_index=True)

    df_gene_exp = df_exp.groupby('gene_id').sum()
    logger.info(f'{df_gene_exp.shape[0]} genes, {df_gene_exp.shape[1]} samples retained for eQTL')

    multi_iso_genes = df_exp['gene_id'].value_counts()
    multi_iso_genes = multi_iso_genes[multi_iso_genes > 1].index
    df_iso_exp = df_exp[df_exp['gene_id'].isin(multi_iso_genes)]
    logger.info(f'{df_iso_exp.shape[0]} isoforms retained (multi-isoform genes only) for isoQTL')

    return df_gene_exp, df_iso_exp


class CallNorm(object):
    def __init__(self, exp, cov, isratio):
        self.df_exp = exp
        self.df_cov = cov
        self.isratio = isratio
        self.df_pheo = None
        self.main()

    def zscore(self, x):
        x = pd.Series(x)
        return stats.norm.ppf((x.rank() - 0.5) / (~pd.isnull(x)).sum())

    def norm(self, df):
        df_in = df.copy().iloc[:, 1:]
        df_in.fillna(0, inplace=True)
        df_in.loc[:, :] = [self.zscore(row) for row in df_in.values]
        return df_in.T

    def splice_ratio(self):
        df = self.df_exp
        ratio = df.iloc[:, 1:] / df.groupby('gene_id')[df.columns[1:]].transform('sum')
        return self.norm(ratio)

    # def lm_covariates(self, exp_mtx, covariates, phenos=None, covs=None):
    #     phenos = phenos or list(exp_mtx.columns)
    #     covs = covs or list(covariates.columns)

    #     logger.info(f'Expression samples: {exp_mtx.shape[0]} | Covariate samples: {covariates.shape[0]}')

    #     exp_mtx = exp_mtx.merge(covariates, left_index=True, right_index=True)
    #     #logger.info(f'Merged samples: {exp_mtx.shape[0]}')
    #     logger.info(f'Starting pre adjust ...')
        
        
    #     resid_infos = [
    #         sm.OLS(exp_mtx[pheno], sm.add_constant(exp_mtx[covs])).fit().resid
    #         for pheno in phenos
    #     ]
    #     return pd.DataFrame(resid_infos, index=phenos, columns=exp_mtx.index).T

    def lm_covariates(self, exp_mtx, covariates, phenos=None, covs=None, block_size=5000):
        phenos = phenos or list(exp_mtx.columns)
        covs = covs or list(covariates.columns)

        logger.info(f'Expression samples: {exp_mtx.shape[0]} | Covariate samples: {covariates.shape[0]}')
        
        # 合并表达矩阵和协变量
        exp_mtx = exp_mtx.merge(covariates, left_index=True, right_index=True)
        logger.info(f'Starting pre-adjust with {len(phenos)} phenotypes, processing in blocks of {block_size} ...')

        # 准备协变量矩阵 X，加截距项
        X = sm.add_constant(exp_mtx[covs].values)  # shape: (n_samples, n_covs + 1)

        # 初始化列表存放所有残差
        resid_infos = []

        # 分块处理 phenos
        for i in range(0, len(phenos), block_size):
            block_phenos = phenos[i:i + block_size]
            logger.info(f'Processing block {i // block_size + 1}: phenos {i} to {min(i + block_size, len(phenos))}')

            for pheno in block_phenos:
                y = exp_mtx[pheno].values  # shape: (n_samples,)
                beta_hat, *_ = np.linalg.lstsq(X, y, rcond=None)
                resid = y - X @ beta_hat
                resid_infos.append(resid)

        # 构建并返回目标格式的 DataFrame（行是样本，列是表型）
        return pd.DataFrame(resid_infos, index=phenos, columns=exp_mtx.index).T



    def main(self):
        if self.df_exp.shape[1] <= 1:
            raise ValueError("Expression matrix has no sample columns")

        df_norm = self.splice_ratio() if self.isratio else self.norm(self.df_exp)
        self.df_pheo = self.lm_covariates(df_norm, self.df_cov)


def write_file(df, out_path):
    logger.info(f'Writing to: {out_path}')
    df.index = ['-'.join(i.split('-')[0:2]) for i in df.index]
    df.index.name = 'IID'
    df.to_csv(out_path, sep='\t')
    
def write_and_export(norm_result, out_prefix, force=False):
    """
    写出表达矩阵，并导出BOD文件
    """
    out_tsv = f"{out_prefix}.ExpNorm.tsv"
    out_bod = f"{out_prefix}.bod"

    # 写表达矩阵
    if os.path.exists(out_tsv) and not force:
        logger.warning(f"Output file exists: {out_tsv}. Skipping.")
    else:
        write_file(norm_result.df_pheo, out_tsv)

    # 写 .bod 文件
    if os.path.exists(out_bod) and not force:
        logger.warning(f"Output file exists: {out_bod}. Skipping.")
    else:
        logger.warning(f"Starting file format to BOD.")
        exp2BOD(out_tsv, out_prefix)

    

def exp2BOD(efile, outpre):
    osca_bin = str(binfinder.find('./resources/osca'))
    if not common.check_file_exists(
        osca_bin,
        file_description=f"OSCA in :{osca_bin}",
        logger=logger,
        exit_on_error=False
    ):
        osca_bin = download_osca()
    
    cmd = [
        osca_bin,
        '--efile', efile,
        '--gene-expression',
        '--make-bod',
        '--no-fid',
        '--out', outpre
    ]

    logger.info(f"Running OSCA make-bod command: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"OSCA failed with return code {e.returncode}")
        raise RuntimeError("OSCA execution failed") from e

def update_opi(opi_file, df_anno, is_gene=False):
    df_opi = pd.read_csv(opi_file, sep='\t', header=None)
    df_anno_slim = df_anno[['chromsome', 'start', 'gene_id', 'strand']].copy()  # 显式创建副本
    df_anno_slim.loc[:, 'chromsome'] = df_anno_slim['chromsome'].apply(lambda x: x.replace('chr', ''))
    
    if is_gene:
        gene2start = {idx: val['start'].min() for idx, val in df_anno_slim.groupby('gene_id')}
        df_anno_gene = df_anno_slim.reset_index().drop(['transcript_id'], axis=1).drop_duplicates('gene_id')
        df_anno_gene.loc[:, 'start'] = df_anno_gene['gene_id'].apply(lambda x: gene2start[x])
        df_anno_gene.loc[:, 'probe'] = df_anno_gene['gene_id']
        df_opi_new = df_opi.merge(df_anno_gene, left_on=1, right_on='probe', how='left')[['chromsome', 'probe', 'start', 'gene_id', 'strand']]
    else:
        df_anno_slim.loc[:, 'probe'] = df_anno_slim.index
        df_opi_new = df_opi.merge(df_anno_slim, left_on=1, right_on='probe', how='left')[['chromsome', 'probe', 'start', 'gene_id', 'strand']]
    
    os.system(f'cp {opi_file} {opi_file}.bak')
    df_opi_new.to_csv(opi_file, sep='\t', index=False, header=None)


def run_preprocess(isoform, covariates, ref, isoform_ratio, prefix, outdir, tpm_threshold, sample_threshold_ratio, force=False):
    outdir = os.path.abspath(outdir or os.path.dirname(isoform))
    os.makedirs(outdir, exist_ok=True)
    out_bod = os.path.join(outdir, 'BOD_files')
    os.makedirs(out_bod, exist_ok=True)

    logger.info(f'Processing isoform file: {isoform}')

    # Step 1: Check and load input files
    logger.info('Checking and loading input files...')
    df_exp, df_anno, df_cov = check_input_files(isoform, covariates, ref)

    # Step 2: Filter isoform and gene expression data
    logger.info('Filtering isoform and gene expression data based on thresholds...')
    df_gene, df_iso = filtered_isoform(df_exp, df_anno, tpm_threshold, sample_threshold_ratio)

    # Step 3: Isoform abundance normalization
    logger.info('Performing isoform abundance normalization...!!!')
    res_iso_abund = CallNorm(df_iso, df_cov.T, isratio=False)
    out_prefix_iso_abund = os.path.join(out_bod, f'{prefix}.isoform_abundance')
    logger.info(f'Writing isoform abundance output to: {out_prefix_iso_abund}')
    write_and_export(res_iso_abund, out_prefix_iso_abund)
    update_opi(f'{out_prefix_iso_abund}.opi', df_anno)
    logger.info('Finished isoform abundance normalization and formated to BOD file')



    # Step 4: Isoform splice ratio normalization (conditional)
    if isoform_ratio:
        logger.info('Performing isoform splice ratio normalization...')
        res_iso_ratio = CallNorm(df_iso, df_cov.T, isratio=True)
        out_prefix_iso_ratio = os.path.join(out_bod, f'{prefix}.isoform_splice_ratio')
        logger.info(f'Writing isoform splice ratio output to: {out_prefix_iso_ratio}')
        write_and_export(res_iso_ratio, out_prefix_iso_ratio)
        update_opi(f'{out_prefix_iso_ratio}.opi', df_anno)
        logger.info('Finished Isoform splice ratio normalization and formated to BOD file')

    else:
        logger.info('Isoform splice ratio normalization not requested. Skipping.')

    # Step 5: Gene abundance normalization
    logger.info('Performing gene abundance normalization...')
    res_gene = CallNorm(df_gene, df_cov.T, isratio=False)
    out_prefix_gene = os.path.join(out_bod, f'{prefix}.gene_abundance')
    logger.info(f'Writing gene abundance output to: {out_prefix_gene}')
    write_and_export(res_gene, out_prefix_gene)
    update_opi(f'{out_prefix_gene}.opi', df_anno,is_gene=True)
    logger.info('Finished gene abundance normalization and formated to BOD file')

    logger.info('All normalization steps completed.')
         
    
@click.command()
@click.option('--verbose', is_flag=True, help='Enable verbose output')
@click.option('--isoform', '-i', required=True, help='Isoform expression file.')
@click.option('--covariates', required=True, help='Covariate file.')
@click.option('--ref', default='gencode_38', type=click.Choice(['refseq_38', 'gencode_38']), help='Reference database.')
@click.option('--isoform-ratio', is_flag=True, help='Calculate isoform expression splice ratio.')
@click.option('--prefix', default='IsoQ', help='Output file prefix. default: IsoQ')
@click.option('--outdir', default='./workdir', help='Output directory. default: workdir')
@click.option('--tpm-threshold', default=0.1, show_default=True, help='TPM threshold for filtering. default: 0.1')
@click.option('--sample-threshold-ratio', default=0.2, show_default=True, 
              help='Minimum fraction of samples where isoform must pass TPM threshold. default: 0.2')
def preprocess(verbose, isoform, covariates, **kwargs):
    """Preprocess input data for IsoQTL"""

    # 设置日志路径
    log_file = f'{datetime.datetime.now().strftime("%Y-%m-%d")}.isoqtl.preprocess.info.log'

    # 初始化日志（自定义的 setup_logger 里完成 format 和 level 设置）
    common.setup_logger(log_file, verbose)
    
    logger = logging.getLogger(__name__)
    logger.info(f'Project starting\nLog file: {log_file}')
    
    # 调用核心逻辑
    run_preprocess(isoform=isoform, covariates=covariates, **kwargs)

if __name__ == '__main__':
    preprocess()
