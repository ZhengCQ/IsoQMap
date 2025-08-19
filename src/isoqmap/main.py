import click
import importlib.metadata
from isoqmap.commands.download import download
from isoqmap.commands.isoquan import isoquan
from isoqmap.commands.isoqtl import isoqtl

__version__ = importlib.metadata.version("isoqmap")


@click.group()
@click.version_option(__version__, prog_name="IsoQMap")
def cli():
    """Isoform Quantification and QTL mapping"""
    pass

cli.add_command(isoquan)
cli.add_command(download)
cli.add_command(isoqtl)

if __name__ == '__main__':
    cli()
