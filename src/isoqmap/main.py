import click
from isoqmap.commands.isoquan import isoquan
from isoqmap.commands.download import download

@click.group()
def cli():
    """Isoform Quantification and QTL mapping """
    pass

cli.add_command(isoquan)
cli.add_command(download)

if __name__ == '__main__':
    cli()