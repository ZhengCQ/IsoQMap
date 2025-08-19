import click
import logging
import os
import glob

from .preprocess import run_preprocess
from .call import run_osca_task
from .format import run_format
from ...tools import pathfinder, common
from ...tools.downloader import download_reference, download_osca
import configparser


logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)

binfinder = pathfinder.BinPathFinder('isoqmap')


def precheck(ref):
    """
    检查并下载 reference files
    """
    try:
        osca_bin = str(binfinder.find('./resources/osca'))
        logger.info(f"Found OSCA binary at: {osca_bin}")
    except FileNotFoundError:
        logger.error("OSCA binary not found. Please install or download it.")
        download_osca()
        osca_bin = str(binfinder.find('./resources/osca'))
        logger.info(f"Downloaded OSCA binary to: {osca_bin}")

    gene_info_fi = str(binfinder.find(f'./resources/ref/{ref}/transcript_gene_info.tsv.gz'))
    if not common.check_file_exists(
        gene_info_fi,
        file_description=f"Gene annotation file {gene_info_fi}",
        logger=logger,
        exit_on_error=False
    ):
        logger.warning(f"Gene annotation file not found. Downloading for {ref}...")
        download_reference(ref, ['geneinfo'])
        gene_info_fi = str(binfinder.find(f'./resources/ref/{ref}/transcript_gene_info.tsv.gz'))

    gene_bed_fi = str(binfinder.find(f'./resources/ref/{ref}/anno_gene_info.bed'))
    logger.info(f"Reference check completed. Gene BED file: {gene_bed_fi}")

@click.command()
@click.option('--ref', type=click.Choice(['refseq_38', 'gencode_38','pig_110']),
              default='gencode_38', help='Reference transcript')
@click.option('--isoform', '-i', required=True, help='Isoform expression file.')
@click.option('--bfile', required=True, help='Prefix of PLINK binary genotype file.')
@click.option('--covariates', required=True, type=click.Path(exists=True), help='Covariate file.')
@click.option('-c', '--config', type=click.Path(exists=True), help='Configuration file')
@click.option('--outdir', default='./workdir', help='Output directory. default: ./workdir')
@click.option('--outprefix', default='osca_qtl_job', help='Prefix of output file. default: osca_qtl_job')
@click.option('--processes', default=5, type=int,
              help='Number of parallel processes to use, default: 5')
@click.option('--force', is_flag=True, default=False, help='Force overwrite of existing files.')

def pipeline(outdir, ref, bfile, isoform, covariates, config, outprefix, processes, force):
    """
    Run the full IsoQTL pipeline: preprocess → call → format
    """

    logger.info("Starting IsoQTL pipeline...")

    # 检查 reference
    precheck(ref)

    # 第一步 preprocess
    logger.info("[Pipeline] Step 1: Preprocessing data...")

    bod_files = [
        os.path.join(outdir, "BOD_files", f"{outprefix}.gene_abundance.bod"),
        os.path.join(outdir, "BOD_files", f"{outprefix}.isoform_abundance.bod"),
        os.path.join(outdir, "BOD_files", f"{outprefix}.isoform_splice_ratio.bod"),
    ]

    if all(os.path.exists(f) for f in bod_files) and not force:
        logger.info("Preprocess output exists. Skipping preprocessing.")
    else:
        run_preprocess(
            isoform=isoform,
            covariates=covariates,
            ref=ref,
            isoform_ratio=True,
            prefix=outprefix,
            outdir=outdir,
            tpm_threshold=0.1,
            sample_threshold_ratio=0.2,
            force=False
        )
    # 第二步 call
    logger.info("[Pipeline] Step 2: Running IsoQTL calls...")

    call_tasks = [
        {
            "name": "eQTL",
            "befile": os.path.join(outdir, "BOD_files", f"{outprefix}.gene_abundance"),
            "mode": "eqtl",
            "outprefix": f'{outprefix}.gene_abundance.eqtl',
            "pattern": os.path.join(outdir, "QTL_results", f"{outprefix}.gene_abundance.eqtl_10_*.besd"),
        },
        {
            "name": "isoQTL",
            "befile": os.path.join(outdir, "BOD_files", f"{outprefix}.isoform_abundance"),
            "mode": "sqtl",
            "outprefix": f'{outprefix}.isoform_abundance.sqtl',
            "pattern": os.path.join(outdir, "QTL_results", f"{outprefix}.isoform_abundance.sqtl_10_*.besd"),
        },
        {
            "name": "irQTL",
            "befile": os.path.join(outdir, "BOD_files", f"{outprefix}.isoform_splice_ratio"),
            "mode": "sqtl",
            "outprefix": f"{outprefix}.isoform_splice_ratio.sqtl",
            "pattern": os.path.join(outdir, "QTL_results", f"{outprefix}.isoform_splice_ratio.sqtl_10_*.besd"),
        }
    ]

    osca_bin = str(binfinder.find('./resources/osca'))
    for task in call_tasks:
        existing_files = glob.glob(task["pattern"])
        if existing_files and not force:
            logger.info(f"{task['name']} results exist. Skipping call.")
        else:
            run_osca_task(
                osca=osca_bin,
                bfile=bfile,
                befile=task["befile"],
                outdir=os.path.join(outdir, "QTL_results"),
                prefix=task["outprefix"],
                mode=task["mode"],
                config = config,
                ref=ref
            )
    # 第三步 format
    logger.info("[Pipeline] Step 3: Formatting QTL results...")

    # isoQTL & irQTL
    run_format(
        verbose=True,
        infile=os.path.join(outdir, "QTL_results", f"{outprefix}.*.sqtl_10_*_isoform_eQTL_effect.txt"),
        mode="sqtl",
        ref=ref,
        id2rs_file=False,
        id2rs_idname='ID',
        id2rs_rsname='rsid',
        processes=processes
        
    )
    # eQTL
    run_format(
        verbose=True,
        infile=os.path.join(outdir, "QTL_results", f"{outprefix}*eqtl_10_*.besd"),
        mode="eqtl",
        ref=ref,
        id2rs_file=False,
        id2rs_idname='ID',
        id2rs_rsname='rsid',
        processes=processes
        
    )
if __name__ == "__main__":
    pipeline()