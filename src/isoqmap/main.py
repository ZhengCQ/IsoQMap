import click
from .commands.isoquan import isoquan

@click.group()
def cli():
    """Splicing-regulatory Driver Genes Identification Tool"""
    pass

cli.add_command(isoquan)

if __name__ == '__main__':
    cli()