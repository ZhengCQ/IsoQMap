import click
from isoqmap.tools.downloader import download_reference

@click.command()
@click.option('--ref', is_flag=True, help='Download reference files')
@click.option('--version', default='gencode_38', help='Reference version, default: gencode_38 ')
@click.option('--files', default='all', type=str,
              help='Files to download: all, transcript, xdata, or comma-separated list, default: all')
@click.pass_context
def download(ctx, ref, version, files):
    """
    Download reference data with hash checking.
    """
    if not ref:
        click.echo("plese set --ref ")

        click.echo(ctx.get_help())
        ctx.exit()

    files_requested = files.split(',') if files != 'all' else ['all']
    download_reference(version, files_requested)
