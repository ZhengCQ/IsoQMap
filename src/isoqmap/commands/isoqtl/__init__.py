import click
from .preprocess import preprocess
from .call import call
from .format import qtlformat
from .pipeline import pipeline

@click.group()
def isoqtl():
    """IsoQTL-related operations: preprocess, run, format, pipeline."""
    pass

isoqtl.add_command(preprocess)
isoqtl.add_command(call)
isoqtl.add_command(qtlformat)
isoqtl.add_command(pipeline)
