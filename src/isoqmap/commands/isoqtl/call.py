import logging
import subprocess
from pathlib import Path
import os
import configparser
import click

from ...tools import pathfinder, common
from ...tools.downloader import download_osca, download_reference

logger = logging.getLogger(__name__)
binfinder = pathfinder.BinPathFinder('isomap')

def resolve_config(config_path):
    cfg = configparser.ConfigParser()
    if config_path:
        cfg.read(config_path, encoding="utf-8")
    else:
        cfg.read(binfinder.find('./config.ini'), encoding="utf-8")
    return cfg

def run_osca_task(osca, bfile, befile, outdir, prefix, mode, config, bed_file=None):
    assert mode in ['sqtl', 'eqtl'], "mode must be 'sqtl' or 'eqtl'"
    cfg = resolve_config(config)

    befile = Path(befile).resolve()
    outdir = Path(outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    out_prefix = outdir / prefix

    task_num = cfg.getint('osca', 'task_num')
    thread_num = cfg.getint('osca', 'thread_num')

    procs = []
    for task_id in range(1, task_num + 1):
        cmd = [
            osca,
            f"--{mode}",
            "--bfile", str(bfile),
            "--befile", str(befile),
            "--maf", str(cfg.get('osca', 'maf')),
            "--call", str(cfg.get('osca', 'call')),
            "--cis-wind", str(cfg.get('osca', 'cis_wind')),
            "--thread-num", str(thread_num),
            "--task-num", str(task_num),
            "--task-id", str(task_id),
            "--out", str(out_prefix)
        ]
        if mode == "sqtl":
            cmd.append("--to-smr")
            if bed_file:
                cmd += ["--bed", str(bed_file)]

        logger.info(f"[Task {task_id}] Running OSCA: {' '.join(cmd)}")
        procs.append(subprocess.Popen(cmd))

    logger.info("Waiting for all OSCA tasks to finish...")
    for p in procs:
        p.wait()
    logger.info("All OSCA tasks finished.")

def write_script(filename: Path, content: str):
    filename = Path(filename)
    filename.write_text(content)
    print(f"写入脚本: {filename}")

def generate_osca_script(osca, bfile, befile, outdir, prefix, mode, bed_file, config, backend):
    cfg = resolve_config(config)
    task_num = cfg.getint('osca', 'task_num')
    thread_num = cfg.getint('osca', 'thread_num')

    cmd = f"""
osca="{osca}"
bfile="{bfile}"
befile="{befile}"
outdir="{outdir}"
prefix="{prefix}"
task_id=$TASK_ID

$osca --{mode} --bfile $bfile --befile $befile --maf {cfg.get('osca','maf')} --call {cfg.get('osca','call')} \\
--cis-wind {cfg.get('osca','cis_wind')} --thread-num {thread_num} --task-num {task_num} --task-id $task_id"""

    if mode == "sqtl":
        cmd += " --to-smr"
        if bed_file:
            cmd += f" --bed {bed_file}"
    cmd += " --out $outdir/$prefix"

    # SLURM
    slurm_script = f"""#!/bin/bash
#SBATCH --job-name={prefix}
#SBATCH --output={outdir}/{prefix}_%A_%a.out
#SBATCH --error={outdir}/{prefix}_%A_%a.err
#SBATCH --array=1-{task_num}
#SBATCH --cpus-per-task={thread_num}
#SBATCH --mem=8G
#SBATCH --time=12:00:00

{cmd.replace('$TASK_ID', '$SLURM_ARRAY_TASK_ID')}
"""

    # SGE
    sge_script = f"""#!/bin/bash
#$ -N {prefix}
#$ -o {outdir}/{prefix}_$TASK_ID.out
#$ -e {outdir}/{prefix}_$TASK_ID.err
#$ -t 1-{task_num}
#$ -pe smp {thread_num}
#$ -l h_vmem=8G

{cmd.replace('$TASK_ID', '$SGE_TASK_ID')}
"""

    # Shell
    shell_script = f"""#!/bin/bash
for task_id in $(seq 1 {task_num}); do
    TASK_ID=$task_id
    {cmd} &
done
wait
"""

    outdir = Path(outdir)
    if backend == 'slurm':
        write_script(outdir / f"run_{prefix}.slurm", slurm_script)
    elif backend == 'sge':
        write_script(outdir / f"run_{prefix}.sge", sge_script)
    elif backend == 'shell':
        write_script(outdir / f"run_{prefix}.sh", shell_script)
    else:
        raise ValueError("backend 必须是 'slurm', 'sge', 或 'shell'")

def batch_generate_scripts(osca, bfile, befile, outdir, prefix, mode, bed_file, config, backend):
    generate_osca_script(osca, bfile, befile, outdir, prefix, mode, bed_file, config, backend)

@click.command()
@click.option('--osca', default=None, help='Path to OSCA binary')
@click.option('--bfile', required=True, help='Prefix for SNP PLINK bfile')
@click.option('--befile', required=True, help='BOD file (expression)')
@click.option('--mode', required=True, type=click.Choice(['sqtl', 'eqtl']), help='QTL analysis mode')
@click.option('--ref', default='gencode_38', type=click.Choice(['refseq_38', 'gencode_38']), help='Reference database')
@click.option('-c', '--config', type=click.Path(exists=True), help='Configuration file')
@click.option('--outdir', default='./workdir', help='Output directory')
@click.option('--prefix', default='osca_qtl_job', help='Output file prefix')
@click.option('--backend', default='shell', type=click.Choice(['slurm', 'sge', 'shell']), help='Execution backend')
@click.option('--run', is_flag=True, help='Whether to run directly')
def call(osca, bfile, befile, mode, ref, outdir, prefix, backend, config, run):
    """Run OSCA QTL Analysis and Generate Job Scripts"""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if osca is None:
        osca = str(binfinder.find('./resources/osca'))
        if not common.check_file_exists(osca, f"OSCA in :{osca}", logger, exit_on_error=False):
            osca = download_osca()

    gene_bed_fi = str(binfinder.find(f'./resources/ref/{ref}/anno_gene_info.bed'))
    if not Path(gene_bed_fi).exists():
        gene_info_fi = str(binfinder.find(f'./resources/ref/{ref}/transcript_gene_info.tsv.gz'))
        if not common.check_file_exists(gene_info_fi, f"Gene annotation file {gene_info_fi}", logger, exit_on_error=False):
            print(f"Gene annotation file not found. Downloading for {ref}...")
            download_reference(ref, ['geneinfo'])
            gene_info_fi = str(binfinder.find(f'./resources/ref/{ref}/transcript_gene_info.tsv.gz'))
        gene_bed_fi = common.geneinfo_2bed(gene_info_fi)
        if not Path(gene_bed_fi).exists():
            raise FileNotFoundError("BED file generation failed.")
    befile_names = os.path.basename(befile)
    prefix = f'{prefix}.{befile_names}'
    bed_file = gene_bed_fi if mode == 'sqtl' else None
    if not prefix.endswith(mode):
        prefix = f"{prefix}.{mode}"

    if run:
        run_osca_task(osca, bfile, befile, outdir, prefix, mode, config, bed_file)
    batch_generate_scripts(osca, bfile, befile, outdir, prefix, mode, bed_file, config, backend)

