import click
from isoqmap.commands.download import download
from isoqmap.commands.isoquan import isoquan
from isoqmap.commands.isoqtl import isoqtl

@click.group()
def cli():
    """Isoform Quantification and QTL mapping"""
    pass

cli.add_command(isoquan)
cli.add_command(download)
cli.add_command(isoqtl)

if __name__ == '__main__':
    cli()

