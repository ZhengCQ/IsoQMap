import click
from isoqmap.commands.isoquan import isoquan
from isoqmap.commands.download import download
from isoqmap.commands.isoqtl_preprocess import preprocess as qtlpreprocess

@click.group()
def cli():
    """Isoform Quantification and QTL mapping"""
    pass  # 不在这里全局检查依赖

cli.add_command(isoquan)
cli.add_command(download)
cli.add_command(qtlpreprocess)

if __name__ == '__main__':
    cli()
