import click
from .preprocess import run_preprocess
from .call import run_osca_task
from .format import run_format
from ...tools import pathfinder,common
import logging

logger = logging.getLogger(__name__)
binfinder = pathfinder.BinPathFinder('isomap')
from ...tools.downloader import download_reference, download_osca


def precheck(ref):
    osca_bin = str(binfinder.find('./resources/osca'))
    
    gene_info_fi = str(binfinder.find(f'./resources/ref/{ref}/transcript_gene_info.tsv.gz')) 
    if not common.check_file_exists(
        gene_info_fi,
        file_description=f"Gene annotaion file {gene_info_fi}",
        logger=logger,
        exit_on_error=False
    ):
        print(f"Gene annotaion file not found or unreadable. Trying to download for {ref}...")
        download_reference(ref, ['geneinfo'])
    
        gene_info_fi = str(binfinder.find(f'./resources/ref/{ref}/transcript_gene_info.tsv.gz'))
           
    gene_bed_fi = str(binfinder.find(f'./resources/ref/{ref}/anno_gene_info.bed'))
    


@click.command()

def pipeline():
    """Run the full IsoQTL pipeline: preprocess -> run -> format"""
    click.echo("Running full IsoQTL pipeline...")

    # 模拟串行调用流程
    click.echo("[Pipeline] Step 1: Preprocessing...")
    run_preprocess()

    click.echo("[Pipeline] Step 2: Running IsoQTL...")
    run_osca_task()

    click.echo("[Pipeline] Step 3: Formatting results...")
    run_format()

    click.echo("Pipeline completed.")



