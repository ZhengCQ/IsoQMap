import os
import subprocess
import logging
from pathlib import Path
import click
from ..tools import pathfinder,common

binfinder = pathfinder.BinPathFinder('isomap')


def run_osca_task(osca_bin, bfile, befile, outdir, prefix, mode, bed=None, task_num=10, threads=4, force=False):
    assert mode in ['sqtl', 'eqtl'], "mode must be 'sqtl' or 'eqtl'"

    befile = Path(befile).resolve()
    outdir = Path(outdir).resolve()
    os.makedirs(outdir, exist_ok=True)
    out_prefix = outdir / prefix

    for task_id in range(1, task_num + 1):
        output_file = f"{out_prefix}.ciseQTL.chr1.txt"  # 可根据实际输出修改
        if not force and Path(output_file).exists():
            logging.info(f"[Task {task_id}] Output exists: {output_file}, skipping.")
            continue

        cmd = [
            osca_bin,
            f"--{mode}",
            "--bfile", str(bfile),
            "--befile", str(befile),
            "--maf", "0.01",
            "--call", "0.85",
            "--cis-wind", "1000",
            "--thread-num", str(threads),
            "--task-num", str(task_num),
            "--task-id", str(task_id),
            "--out", str(out_prefix)
        ]
        if mode == "sqtl":
            cmd.append("--to-smr")
            if bed:
                cmd += ["--bed", bed]

        logging.info(f"[Task {task_id}] Running OSCA: {' '.join(cmd)}")
        subprocess.Popen(cmd)

    logging.info("Waiting for all OSCA tasks to finish...")
    os.wait()
    logging.info("All OSCA tasks finished.")

def write_script(filename, content):
    with open(filename, 'w') as f:
        f.write(content)
    print(f"✅ 写入脚本: {filename}")

def generate_osca_script(prefix, mode, befile, outdir, bedfile=None, osca_path='osca',
                         bfile='data/WGS/bedfile/snp.filtered.noPar',
                         task_num=10, threads=4, backend='slurm'):
    cmd = f"""osca="{osca_path}"
bfile="{bfile}"
befile="{befile}"
outdir="{outdir}"
prefix="{prefix}"
task_id=$TASK_ID

$osca --{mode} --bfile $bfile --befile $befile --maf 0.01 --call 0.85 \\
--cis-wind 1000 --thread-num {threads} --task-num {task_num} --task-id $task_id"""

    if mode == "sqtl":
        cmd += " --to-smr"
        if bedfile:
            cmd += f" --bed {bedfile}"
    cmd += " --out $outdir/$prefix"

    # 各种脚本模板
    slurm_script = f"""#!/bin/bash
#SBATCH --job-name={prefix}
#SBATCH --output={outdir}/{prefix}_%A_%a.out
#SBATCH --error={outdir}/{prefix}_%A_%a.err
#SBATCH --array=1-{task_num}
#SBATCH --cpus-per-task={threads}
#SBATCH --mem=8G
#SBATCH --time=12:00:00

{cmd.replace('$TASK_ID', '$SLURM_ARRAY_TASK_ID')}
"""

    sge_script = f"""#!/bin/bash
#$ -N {prefix}
#$ -o {outdir}/{prefix}_$TASK_ID.out
#$ -e {outdir}/{prefix}_$TASK_ID.err
#$ -t 1-{task_num}
#$ -pe smp {threads}
#$ -l h_vmem=8G

{cmd.replace('$TASK_ID', '$SGE_TASK_ID')}
"""

    shell_script = f"""#!/bin/bash
# 普通 shell 执行脚本，循环运行所有任务
for task_id in $(seq 1 {task_num}); do
    TASK_ID=$task_id
    {cmd}
done
"""

    if backend == 'slurm':
        write_script(f"run_{prefix}.slurm", slurm_script)
    elif backend == 'sge':
        write_script(f"run_{prefix}.sge", sge_script)
    elif backend == 'shell':
        write_script(f"run_{prefix}.sh", shell_script)
    else:
        raise ValueError("backend 必须是 'slurm', 'sge', 或 'shell'")

def batch_generate_scripts(prefix, mode, befile, outdir, bed_file, osca_path, bfile, backend):
    generate_osca_script(prefix, mode, befile, outdir, bed_file,
                         osca_path=osca_path, bfile=bfile, backend=backend)

@click.command()
@click.option('--osca', default=None, help='Path to OSCA binary')
@click.option('--bfile', required=True, help='Prefix for SNP PLINK bfile')
@click.option('--befile', required=True, help='BOD file (expression)')
@click.option('--mode', required=True, type=click.Choice(['sqtl', 'eqtl']), help='QTL analysis mode')
@click.option('--outdir', default='osca_out', help='Output directory')
@click.option('--prefix', default='osca_job', help='Output file prefix')
@click.option('--backend', default='shell', type=click.Choice(['slurm', 'sge', 'shell']), help='Execution backend')
@click.option('--run', is_flag=True, help='Whether to run directly in Python')
def main(osca, bfile, befile, mode, outdir, prefix, backend, run):
    from pathlib import Path

    # 自动寻找 OSCA 二进制文件
    if osca is None:
        osca = binfinder('./resources/osca')
        if osca is None:
            raise FileNotFoundError("未提供 --osca 参数且 ./resources/osca 下未找到 OSCA 可执行文件")

    # bed 文件在 sqtl 模式下才使用
    bed_file = '/share/Apps/iGTEx/ref/gencode_38/anno_gene_info.bed' if mode == 'sqtl' else None

    if run:
        run_osca_task(osca, bfile, befile, outdir, prefix, mode, bed_file)
    else:
        batch_generate_scripts(prefix, mode, befile, outdir, bed_file, osca_path=osca, bfile=bfile, backend=backend)

