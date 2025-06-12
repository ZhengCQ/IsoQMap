import click
from .preprocess import run_preprocess
from .call import run_osca_task
from .format import run_format

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
