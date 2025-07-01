import click
import logging
import pandas as pd
from pathlib import Path
import glob
from multiprocessing import Pool
import numpy as np
import gzip
import datetime
import subprocess
import sys
import os

from ...tools import pathfinder, common
from ...tools.downloader import download_reference,download_osca

# Initialize path finder
binfinder = pathfinder.BinPathFinder('isomap')

# Global variables
p_cut = 5e-8
id2rs = {}
iso2tss = {}
gene2tss = {}
logger = logging.getLogger(__name__)  # Define logger at module level

def load_global_data(anno_path, geneinfo_path, id_col='ID', rsid_col='rsid'):
    """Load global data required for processing"""
    global id2rs, iso2tss, gene2tss
    
    try:
        # Load SNP annotation data
        if anno_path:
            logger.info(f"Loading SNP annotation data from {anno_path}")
            df_anno = pd.read_csv(anno_path, sep='\t')
            df_anno = df_anno[~pd.isnull(df_anno[rsid_col])]
            id2rs = {i[0]: i[1] for i in df_anno[[id_col, rsid_col]].to_dict('split')['data']}
            logger.info(f"Loaded {len(id2rs)} SNP mappings")
        else:
            logger.info(f"Not exists anno_path, ingore id2rs")

        # Load gene and transcript information
        logger.info(f"Loading gene info from {geneinfo_path}")
        df_geneinfo = pd.read_csv(geneinfo_path, sep='\t')
        iso2tss = {i[1]: int(i[0]) for i in df_geneinfo[['start', 'transcript_id']].to_dict('split')['data']}
        gene2tss = {idx: int(val['start'].min()) for idx, val in df_geneinfo[['start', 'gene_id']].groupby('gene_id')}
        logger.info(f"Loaded {len(gene2tss)} gene mappings and {len(iso2tss)} isoform mappings")
    
    except Exception as e:
        logger.error(f"Error loading global data: {str(e)}")
        raise

def format_file(fi):
    """Process input file and generate formatted output files"""
    header = ['SNP', 'chr', 'pos', 'A1', 'A2', 'freq', 'beta', 'se', 'pval', 'pheno_name', 'tss_dis']
    try:
        if not fi.endswith('gz'):
            base_name = fi.replace(".txt", "")
            cmd = f'gzip -f {fi}'
            logger.info(f"Running gzip command: {cmd}")
            #subprocess.run(cmd, check=True)
            os.system(f'{cmd}')
            fi = f'{fi}.gz'
        else:
            base_name = fi.replace(".txt.gz", "")
            
        logger.info(f"Formatting file: {fi}")
        
        # Open output files
        with open(f"{base_name}.format.gene.txt", 'w') as outf_gene, \
             open(f"{base_name}.format.isoform.txt", 'w') as outf_iso:
            
            # Write headers to all output files
            outf_gene.write('\t'.join(header) + '\n')
            outf_iso.write('\t'.join(header + ['pheno_genes']) + '\n')

            # Process input file
            with gzip.open(fi, 'rt') as f:
                next(f)  # Skip header
                for line in f:
                    items = line.strip().split('\t')
                    
                    # Extract common fields
                    ids = items[0]
                    chrom = items[1]
                    pos = items[2]
                    a1 = items[3]
                    a2 = items[4]
                    freq = items[5]
                    beta = items[10]
                    se = items[11]
                    p = items[12]
                    pheno_name_gene = items[6]
                    
                    if p is None:
                        continue

                    # Apply filters
                    if float(freq) < 0.05:
                        continue
                    
                    if not id2rs:
                        snp = ids
                    else:
                        if ids not in id2rs:
                            snp = ids
                        else:
                            snp = id2rs[ids]                        
                    
                    if pheno_name_gene not in gene2tss:
                        continue
                    
                    tss_dis = gene2tss[pheno_name_gene] - int(pos)
                    
                    if abs(tss_dis) > 1000000:
                        continue
                    
                    # Gene-level output
                    gene_output = '\t'.join([
                        snp, chrom, pos, a1, a2, freq, beta, se, p, 
                        pheno_name_gene, str(tss_dis)
                    ]) + '\n'
                    outf_gene.write(gene_output)
                    
                    # Isoform-level processing
                    isf_eqtl = np.array(items[13:]).reshape(-1, 4)
                    for each in isf_eqtl:
                        pheno_name_t = each[0]
                        beta_t = each[1]
                        se_t = each[2]
                        p_t = each[3]
                        
                        if pheno_name_t not in iso2tss:
                            continue
                        
                        tss_dis_t = iso2tss[pheno_name_t] - int(pos)
                        isoform_output = '\t'.join([
                            snp, chrom, pos, a1, a2, freq, 
                            beta_t, se_t, p_t, pheno_name_t, str(tss_dis_t), pheno_name_gene
                        ]) + '\n'
                        outf_iso.write(isoform_output)

    except Exception as e:
        logger.error(f"Error processing file {fi}: {str(e)}")
        raise

def fetch_sig(fi):
    """Extract significant results from formatted files"""
    try:
        if fi.endswith('.txt.gz'):
            base_name = fi.replace(".txt.gz", "")
        else:
            base_name = fi.replace(".txt", "")
        logger.info(f"Extracting significant results from {fi}")
        
        # For Gene
        fi_gene = f'{base_name}.format.gene.txt'
        df_format = pd.read_csv(fi_gene, sep='\t')
        df_format = df_format.dropna()
        
        # Get significant genes (p < cutoff)
        sig_genes = df_format[df_format['pval'] < p_cut]['pheno_name'].drop_duplicates()
        
        # Write various output files
        outputs = [
            (df_format[df_format['pheno_name'].isin(sig_genes)], 
             f'{base_name}.format.sigGene.txt.gz'),
             
            (df_format.loc[df_format.groupby('pheno_name')['pval'].idxmin()], 
             f'{base_name}.format.leadSNP_gene.txt.gz'),
             
            (df_format[df_format['pval'] < p_cut], 
             f'{base_name}.format.gene.signifpairs.txt.gz')
        ]
        
        for df, outfile in outputs:
            df.to_csv(outfile, sep='\t', index=False)
            logger.info(f"Saved {outfile}")

        # For isoform
        fi_isoform = f'{base_name}.format.isoform.txt'
        df_format_isoform = pd.read_csv(fi_isoform, sep='\t')
        df_format_isoform = df_format_isoform.dropna()

        sig_isoform = df_format_isoform[df_format_isoform['pval'] < p_cut]['pheno_name'].drop_duplicates()
        df_format_isoform_sig = df_format_isoform[
            (df_format_isoform['pheno_genes'].isin(sig_genes)) & 
            (df_format_isoform['pheno_name'].isin(sig_isoform))
        ]
        
        iso_outputs = [
            (df_format_isoform_sig, f'{base_name}.format.sigIsoform.txt.gz'),
            (df_format_isoform.loc[df_format_isoform.groupby(['pheno_genes','pheno_name'])['pval'].idxmin()], 
             f'{base_name}.format.leadSNP_isoform.txt.gz'),
            (df_format_isoform[df_format_isoform['pval'] < p_cut], 
             f'{base_name}.format.isoform.signifpairs.txt.gz')
        ]
        
        for df, outfile in iso_outputs:
            df.to_csv(outfile, sep='\t', index=False)
            logger.info(f"Saved {outfile}")
            
        logger.info(f"Compress {fi_gene} and {fi_isoform}")
        os.system(f'gzip {fi_gene}')
        os.system(f'gzip {fi_isoform}')
        logger.info(f"Get {fi_gene}.gz and {fi_isoform}.gz") 
        
    
    except Exception as e:
        logger.error(f"Error processing significant results for {fi}: {str(e)}")
        raise



def besd2txt(fi):
    """Convert BESD format to text format using OSCA"""
    try:
        osca_bin = str(binfinder.find('./resources/osca'))
        if not common.check_file_exists(
            osca_bin,
            file_description=f"OSCA in :{osca_bin}",
            logger=logger,
            exit_on_error=False
        ):
            osca_bin = download_osca()

        basename = fi.replace(".besd", "")
        output_file = f'{basename}.txt'

        if os.path.exists(f"{output_file}.gz"):
            os.remove(f"{output_file}.gz")

        # 运行 OSCA
        cmd_osca = [osca_bin, "--beqtl-summary", basename, "--query", "1", "--out", output_file]
        logger.info(f"Running OSCA: {' '.join(cmd_osca)}")
        subprocess.run(cmd_osca, check=True)

        # 确认输出文件存在后再压缩
        if os.path.exists(output_file):
            cmd_gzip = ["gzip", "-f", output_file]
            logger.info(f"Compressing: {' '.join(cmd_gzip)}")
            subprocess.run(cmd_gzip, check=True)
        else:
            raise FileNotFoundError(f"{output_file} not generated by OSCA!")
    except subprocess.CalledProcessError as e:
        logger.error(f"OSCA failed with return code {e.returncode}")
        raise RuntimeError("OSCA execution failed") from e
    except Exception as e:
        logger.error(f"Error in besd2txt for file {fi}: {str(e)}")
        raise

def process_gene_data(input_path: str) -> None:
    """Process single gene file (multi-processing compatible version)"""
    try:
        
        output_path = Path(input_path).with_name(
            Path(input_path).name.replace(".txt.gz", ".format.tsv.gz")
        )
        
        logger.info(f"Processing: {input_path} → {output_path}")
        
        # Read data
        df_gene = pd.read_csv(input_path, sep='\t', compression='gzip')
        
        # Validate required columns
        required_columns = {'SNP', 'Gene', 'BP', 'Freq', 'b', 'SE', 'p', 'Probe'}
        if not required_columns.issubset(df_gene.columns):
            missing = required_columns - set(df_gene.columns)
            raise ValueError(f"Missing columns: {missing}")

        # Process data
        if id2rs:
            df_gene['SNP'] = df_gene['SNP'].apply(lambda x:id2rs.get(x, x))
            
        df_gene = df_gene[df_gene['Gene'].isin(gene2tss)]
        df_gene['tss_dis'] = df_gene['Gene'].map(gene2tss) - df_gene['BP']
        
        df_gene = df_gene[(df_gene['tss_dis'].abs() <= 1_000_000) & (df_gene['Freq'] > 0.05)]
        
        # Rename and select columns
        df_gene_slim = df_gene.rename(columns={
            'Chr': 'chr', 'BP': 'pos', 'b': 'beta', 'Freq': 'freq',
            'SE': 'se', 'p': 'pval', 'Probe': 'pheno_name'
        })[['SNP', 'chr', 'pos', 'A1', 'A2', 'freq', 'beta', 'se', 'pval', 'pheno_name', 'tss_dis']]
        
        # Save results
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df_gene_slim.to_csv(output_path, sep='\t', index=False)
        
        # Generate additional output files
        sig_genes = df_gene_slim[df_gene_slim['pval'] < 5e-8]['pheno_name'].drop_duplicates()
        
        outputs = [
            (df_gene_slim[df_gene_slim['pheno_name'].isin(sig_genes)],
             str(output_path).replace('.format.tsv.gz', '.format.sigGene.tsv.gz')),
             
            (df_gene_slim.loc[df_gene_slim.groupby('pheno_name')['pval'].idxmin()],
             str(output_path).replace('.format.tsv.gz', '.format.leadSNP_gene.tsv.gz')),
             
            (df_gene_slim[df_gene_slim['pval'] < 5e-8],
             str(output_path).replace('.format.tsv.gz', '.format.signifpairs.tsv.gz'))
        ]

        for df, out_path in outputs:
            df.to_csv(out_path, sep='\t', index=False)
            logger.info(f"Saved {out_path}")

    except Exception as e:
        logger.error(f"Failed to process {input_path}: {str(e)}")
        raise

def safe_format_file(fi):
    """Wrapper function to handle exceptions in format_file"""
    try:
        format_file(fi)
    except Exception as e:
        logger.error(f"Failed to process file {fi}: {str(e)}")
        return False
    return True

def safe_fetch_sig(fi):
    """Wrapper function to handle exceptions in fetch_sig"""
    try:
        fetch_sig(fi)
    except Exception as e:
        logger.error(f"Failed to extract significant results from {fi}: {str(e)}")
        return False
    return True

def safe_besd2txt(fi):
    """Wrapper function to handle exceptions in besd2txt"""
    try:
        besd2txt(fi)
    except Exception as e:
        logger.error(f"Failed to convert BESD file {fi}: {str(e)}")
        return False
    return True

def safe_process_gene_data(fi):
    """Wrapper function to handle exceptions in process_gene_data"""
    try:
        process_gene_data(fi)
    except Exception as e:
        logger.error(f"Failed to process gene data file {fi}: {str(e)}")
        return False
    return True

def run_format(verbose, infile, mode, ref, id2rs_file, id2rs_idname, id2rs_rsname, processes):
    """Main processing function"""
    try:
        # Find gene info file
        gene_info_fi = str(binfinder.find(f'./resources/ref/{ref}/transcript_gene_info.tsv.gz'))    
        
        if not common.check_file_exists(
            gene_info_fi,
            file_description=f"Gene annotation file {gene_info_fi}",
            logger=logger,
            exit_on_error=False
        ):
            logger.info(f"Gene annotation file not found. Downloading for {ref}...")
            download_reference(ref, ['geneinfo'])
            gene_info_fi = str(binfinder.find(f'./resources/ref/{ref}/transcript_gene_info.tsv.gz')) 

        # Load required data
        load_global_data(id2rs_file, gene_info_fi, id2rs_idname, id2rs_rsname)
        
        # Process files based on mode
        files_to_process = glob.glob(infile)

        
        if not files_to_process:
            logger.error(f"No files matched pattern: {infile}")
            return
        
        logger.info(f"Found {len(files_to_process)} files to process")

        with Pool(processes=processes) as pool:
            if mode == 'sqtl':
                # First pass: format files
                results = pool.map(safe_format_file, files_to_process)
                failed = sum(1 for r in results if not r)
                if failed > 0:
                    logger.warning(f"Failed to process {failed} files in format step")
                
                # Second pass: extract significant results
                results = pool.map(safe_fetch_sig, files_to_process)
                failed = sum(1 for r in results if not r)
                if failed > 0:
                    logger.warning(f"Failed to extract significant results from {failed} files")
            else:
                ## check and download osca
                osca_bin = str(binfinder.find('./resources/osca'))
                if not common.check_file_exists(
                    osca_bin,
                    file_description=f"OSCA in :{osca_bin}",
                    logger=logger,
                    exit_on_error=False
                ):
                    osca_bin = download_osca()
                
                
                files_besd_to_process = []
                for i in files_to_process:
                    if not i.endswith('.besd') and not os.path.exists(i):
                        logger.warning(f"{i} failed to query from besd to txt because {i} not exits or not endwiths *besd ")
                        continue
                    elif os.path.exists(i.replace('.besd', '.txt.gz')):
                        logger.warning(f"{i} failed to query from besd to txt because results {i.replace('.besd', '.txt.gz')} exits")
                        continue
                    files_besd_to_process.append(i)  
                
                if len(files_besd_to_process)>0:
                    # convert BESD files
                    results = pool.map(safe_besd2txt, files_besd_to_process)
                    failed = sum(1 for r in results if not r)
                    if failed > 0:
                        logger.warning(f"Failed to convert {failed} BESD files")
                              
                # Second pass: process gene data
                txt_files = []  
                for i in files_to_process:
                    txt_fi = i.replace('.besd', '.txt.gz')
                    if not os.path.exists(txt_fi) or not txt_fi.endswith('.txt.gz'):
                        logger.warning(f"{txt_fi} failed to query from besd to txt because {txt_fi} not exits or not endwiths *besd ")
                        continue
                    txt_files.append(txt_fi)
                results = pool.map(safe_process_gene_data, txt_files)
                failed = sum(1 for r in results if not r)
                if failed > 0:
                    logger.warning(f"Failed to process {failed} gene data files")
        
        logger.info("Processing completed with %d errors" % failed)

    except Exception as e:
        logger.error(f"Error in run_format: {str(e)}")
        raise

@click.command(name="format")
@click.option('--verbose', is_flag=True, help='Enable verbose output')
@click.option('--infile', required=True, help='Input osca output file. For eqtl: *besd file (e.g. gene_abudance.eqtl_10_*.besd); For sqtl, *isoform_eQTL_effect.txt (e.g. isoform_splice_ratio.sqtl_10-*.isoform_eQTL_effect.txt)')
@click.option('--mode', required=True, type=click.Choice(['sqtl', 'eqtl']), 
              help='QTL analysis mode (sQTL or eQTL)')
@click.option('--id2rs-file', required=False, 
              help='Path to file mapping variant IDs to rsIDs')
@click.option('--ref', default='gencode_38', 
              type=click.Choice(['refseq_38', 'gencode_38']), 
              help='Reference database')
@click.option('--id2rs-idname', default='ID', 
              help='Column name for variant ID in id2rs file, defualt: ID')
@click.option('--id2rs-rsname', default='rsid', 
              help='Column name for rsID in id2rs file, defualt: rsid')
@click.option('--processes', default=10, type=int,
              help='Number of parallel processes to use, default: 10')

def qtlformat(verbose, infile, mode, **kwargs):
    """Format QTL results for downstream usage"""
    # Set up logging
    log_file = f'{datetime.datetime.now().strftime("%Y-%m-%d")}.isoqtl.format.info.log'

    common.setup_logger(log_file, verbose)
    
    global logger
    logger = logging.getLogger(__name__)
    logger.info(f'Project starting\nLog file: {log_file}')
    
    # Run main processing
    try:
        run_format(infile=infile, mode=mode, verbose=verbose, **kwargs)
    except Exception as e:
        logger.error(f"Fatal error: {str(e)}")
        sys.exit(1)

if __name__ == '__main__':
    qtlformat()