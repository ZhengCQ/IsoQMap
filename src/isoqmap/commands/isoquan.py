#!/usr/bin/env python

import os
import threading
import logging
import datetime
import configparser
import subprocess
import concurrent.futures
import pandas as pd
import click
import queue
from ..tools import pathfinder, common
from ..tools.downloader import download_reference

binfinder = pathfinder.BinPathFinder('isoqmap')

# 全局线程锁
threadLock = threading.Lock()

# ----------------------------------------------------------------------------
# Logging utils
# ----------------------------------------------------------------------------

def setup_logger(log_file, verbose=False):
    FORMAT = '%(asctime)s %(message)s'
    level = logging.DEBUG if verbose else logging.INFO

    logger = logging.getLogger()

    if logger.hasHandlers():
        logger.handlers.clear()

    logging.basicConfig(filename=log_file, level=level, format=FORMAT)

    console = logging.StreamHandler()
    console.setLevel(level)
    formatter = logging.Formatter(FORMAT)
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)


def log_info(msg):
    logging.getLogger(__name__).info(msg)

def log_error(msg):
    logging.getLogger(__name__).error(msg)

# ----------------------------------------------------------------------------
# File utilities
# ----------------------------------------------------------------------------

def write_shell(shell_file, cmd):
    os.makedirs(os.path.dirname(shell_file), exist_ok=True)
    with open(shell_file, 'w') as f:
        f.write(cmd)

# ----------------------------------------------------------------------------
# Class for job status
# ----------------------------------------------------------------------------

class JobStatus:
    def __init__(self, sh_nm, df_status):
        self.sh_nm = sh_nm
        self.df_status = df_status
        self.run()

    def run(self):
        log_info(f"Running: {self.sh_nm}")
        res = subprocess.run(
            f'bash {self.sh_nm} 1>{self.sh_nm}.stdout 2>{self.sh_nm}.stderr',
            shell=True
        )
        if res.returncode == 0:
            log_info(f"Success: {self.sh_nm}")
            self.change_status('Success')
        else:
            log_error(f"Error: Please check {self.sh_nm}.stderr")
            self.change_status('Error')

    def change_status(self, status_item):
        with threadLock:
            self.df_status.loc[self.df_status['shell'] == self.sh_nm, 'status'] = status_item

# ----------------------------------------------------------------------------
# Thread class
# ----------------------------------------------------------------------------

class MyThread(threading.Thread):
    def __init__(self, q, df_status):
        threading.Thread.__init__(self)
        self.q = q
        self.df_status = df_status

    def run(self):
        while True:
            try:
                cmd = self.q.get(timeout=2)
                job = JobStatus(cmd, self.df_status)
                self.df_status = job.df_status
            except queue.Empty:
                break

# ----------------------------------------------------------------------------
# Functional blocks
# ----------------------------------------------------------------------------

def check_fq(fqfile):
    suffixs = ('fq.gz', 'fq', 'fastq', 'fastq.gz', 'fa.gz', 'fa', 'fasta', 'fasta.gz')
    if not fqfile.endswith(suffixs):
        log_error(f"This is not an fq/fa file. Suffix should be one of: {suffixs}")
        raise click.BadParameter(f"Invalid file format: {fqfile}")

def read_sampleinfo(infile):
    sample_info = []
    with open(infile) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            sample, lib, fq1, fq2 = parts[:4]
            for fq in (fq1, fq2):
                if not os.path.exists(fq):
                    log_error(f"{fq} in {sample} not exists!")
                    raise click.BadParameter(f"File not found: {fq}")
            sample_info.append([sample, lib, os.path.abspath(fq1), os.path.abspath(fq2)])
    return pd.DataFrame(sample_info, columns=['sample', 'lib', 'fq1', 'fq2'])

def ensure_transcript_exists(refdb, config, binfinder):
    transcript = config.get('xaem', 'transcript_fa') or str(
        binfinder.find(f'./resources/ref/{refdb}/transcript.fa.gz')
    )
    if not common.check_file_exists(
        transcript,
        file_description=f"transcript file is {transcript} for {refdb}",
        logger=logging.getLogger(__name__),
        exit_on_error=False
    ):
        log_info(f"Transcript file not found. Trying to download for {refdb}...")
        download_reference(refdb, ['transcript'])
        transcript = str(binfinder.find(f'./resources/ref/{refdb}/transcript.fa.gz'))
        if not common.check_file_exists(
            transcript,
            file_description=f"transcript file is {transcript} for {refdb}",
            logger=logging.getLogger(__name__),
            exit_on_error=True
        ):
            raise FileNotFoundError(f"Transcript still not found after download for {refdb}")
    return transcript

def index_ref(outdir, config, xaem_dir, refdb, step=1):
    transcript = ensure_transcript_exists(refdb, config, binfinder)
    outfa = f'{outdir}/ref/{os.path.basename(transcript).replace(".gz", "")}'
    index_dir = f'{outdir}/ref/TxIndexer_idx'

    cmd = ""
    if transcript.endswith('gz'):
        cmd += f"gunzip -c {transcript} > {outfa}\n"
    else:
        cmd += f"ln -fs {transcript} {outfa}\n"
    cmd += f"sed -i 's/|/ /' {outfa}\n"
    cmd += f"{xaem_dir}/bin/TxIndexer -t {outfa} --out {index_dir}\n"

    shell_file = f'{outdir}/shell/Step{step}.index_fa.sh'
    write_shell(shell_file, cmd)
    return index_dir, shell_file

def get_eqclass(df, outdir, xaem_dir, TxIndexer_idx, step=2):
    shell_lst = []
    seqdir = f'{outdir}/seqData'
    os.makedirs(seqdir, exist_ok=True)

    for sample, val in df.groupby('sample'):
        fq1_lst = list(val['fq1'])
        fq2_lst = list(val['fq2'])
        sample_fq1 = f'{seqdir}/{sample}_1.fq.gz'
        sample_fq2 = f'{seqdir}/{sample}_2.fq.gz'

        cmd = ""
        if len(fq1_lst) == 1:
            cmd += f"ln -fs {fq1_lst[0]} {sample_fq1}\n"
            cmd += f"ln -fs {fq2_lst[0]} {sample_fq2}\n"
        else:
            cmd += f"zcat {' '.join(fq1_lst)} | gzip -cf > {sample_fq1} &\n"
            cmd += f"zcat {' '.join(fq2_lst)} | gzip -cf > {sample_fq2} &\n"
            cmd += "wait\n"

        cmd += f"""{xaem_dir}/bin/XAEM \\
    -i {TxIndexer_idx} \\
    -l IU \\
    -1 <(gunzip -c {sample_fq1}) \\
    -2 <(gunzip -c {sample_fq2}) \\
    -p 2 \\
    -o {outdir}/results/{sample}\n"""

        shell_file = f'{outdir}/shell/Step{step}.gen_eqclass_{sample}.sh'
        write_shell(shell_file, cmd)
        shell_lst.append(shell_file)
    return shell_lst

def ensure_xmatrix_exists(refdb, config, binfinder):
    xmatrix = config.get('xaem', 'x_matrix') or str(
        binfinder.find(f'./resources/ref/{refdb}/X_matrix.RData')
    )
    if not common.check_file_exists(
        xmatrix,
        file_description=f"X_matrix file is {xmatrix} for {refdb}",
        logger=logging.getLogger(__name__),
        exit_on_error=False
    ):
        log_info(f"X_matrix file not found. Trying to download for {refdb}...")
        download_reference(refdb, ['X_matrix'])
        xmatrix = binfinder.find(f'./resources/ref/{refdb}/X_matrix.RData')
        if not common.check_file_exists(
            xmatrix,
            file_description=f"X_matrix file is {xmatrix} for {refdb}",
            logger=logging.getLogger(__name__),
            exit_on_error=True
        ):
            raise FileNotFoundError(f"X_matrix still not found after download for {refdb}")
    return xmatrix

def count_matrix(outdir, xaem_dir, config, x_matrix, step=3):
    resdir = f'{outdir}/results'
    cmd = "\n".join([
        f"Rscript {xaem_dir}/R/Create_count_matrix.R workdir={resdir} core=8 design.matrix={x_matrix}",
        f"""Rscript {xaem_dir}/R/AEM_update_X_beta.R \\
    workdir={resdir} \\
    core={config.getint('xaem', 'update_cpu')} \\
    design.matrix={x_matrix} \\
    merge.paralogs={config.getboolean('xaem', 'merge.paralogs')} \\
    isoform.method={config.get('xaem', 'isoform.method')} \\
    remove.ycount={config.getboolean('xaem', 'remove.ycount')}""",
        f"Rscript {binfinder.find('./tools/isoform_rdata2exp.R')} inRdata={resdir}/XAEM_isoform_expression.RData && gzip {resdir}/XAEM_isoform_expression_count.tsv && gzip {resdir}/XAEM_isoform_expression_tpm.tsv"
    ])
    shell_file = f'{outdir}/shell/Step{step}.matrix_samples.sh'
    write_shell(shell_file, cmd)
    return shell_file


# -----------------------------------------------------------------------
# Job run helpers
# -----------------------------------------------------------------------

def single_job_run(sh_nm, df, status_file):
    job = JobStatus(sh_nm, df)
    df_status = job.df_status
    write_status(df_status, status_file)
    return df_status

def is_success(df, step_name):
    df_step = df[df['name'] == step_name]
    return df_step.shape[0] == (df_step['status'] == 'Success').sum()

def write_status(df, status_file):
    df.to_csv(status_file, sep='|', index=False)
    df[df['status'] == 'Error'].to_csv(f'{status_file}.Error', sep='|', index=False)


def get_all_shells(outdir, df_sample, config, xaem_dir, refdb, xaem_index=None, x_matrix=None):
    shell_info = []
    step_n = 1
    if xaem_index:
        TxIndexer_idx = os.path.abspath(xaem_index)
    else:
        TxIndexer_idx, cmd = index_ref(outdir, config, xaem_dir, refdb, step=step_n)
        shell_info.append([cmd, step_n, 'index'])
        step_n += 1

    eqclass_shells = get_eqclass(df_sample, outdir, xaem_dir, TxIndexer_idx, step=step_n)
    shell_info.extend([[i, step_n, 'eqclass'] for i in eqclass_shells])
    step_n += 1

    x_matrix = ensure_xmatrix_exists(refdb, config, binfinder)

    cmd = count_matrix(outdir, xaem_dir, config, x_matrix, step=step_n)
    shell_info.append([cmd, step_n, 'matrix'])

    df_shell_info = pd.DataFrame(shell_info)
    df_shell_info['status'] = 'Ready'
    df_shell_info.columns = ['shell', 'step', 'name', 'status']
    status_file = f'{outdir}/shell/JOB.Status'
    df_shell_info.to_csv(status_file, index=False, sep='|')
    return status_file
# -----------------------------------------------------------------------
# Main logic
# -----------------------------------------------------------------------



def run_isoquan(infile, ref, config, outdir, xaem_dir, xaem_index, x_matrix, force):
    """
    Run the full isoform quantification pipeline
    """
    # Load config
    cfg = configparser.ConfigParser()
    if config:
        cfg.read(config, encoding="utf-8")
    else:
        cfg.read(binfinder.find('./config.ini'), encoding="utf-8")

    # XAEM path
    if not xaem_dir:
        xaem_dir = cfg.get('xaem', 'xaem_dir') or binfinder.find('./resources/XAEM/XAEM-binary-0.1.1-cq')
    log_info(f"Parameter: xaem_dir is {xaem_dir}")
    if not os.path.exists(xaem_dir):
        raise FileNotFoundError(f"XAEM binary not found at {xaem_dir}")

    # Prepare directories
    outdir = os.path.abspath(outdir)
    for sub in ["", "seqData", "results", "shell", "ref"]:
        os.makedirs(os.path.join(outdir, sub), exist_ok=True)

    # Sample info
    df_sample = read_sampleinfo(infile)

    # Prepare shell scripts
    status_file = f"{outdir}/shell/JOB.Status"
    if not os.path.exists(status_file) or force:
        log_info("Generating all jobs...")
        status_file = get_all_shells(
            outdir, df_sample, cfg, xaem_dir, ref,
            xaem_index=xaem_index, x_matrix=x_matrix
        )
    else:
        log_info(f"Found existing status file {status_file}. Resuming jobs. Use --force to regenerate.")

    # Load job status
    df_status = pd.read_csv(status_file, sep="|")
    status_dict = df_status.groupby("status").size().to_dict()
    log_info(f"There are {df_status.shape[0]} jobs: {status_dict}")

    # ----------------------
    # Step 1: index
    # ----------------------
    df_index = df_status[df_status['name'] == 'index']
    if not is_success(df_status, "index") and not df_index.empty:
        sh = df_index.iloc[0]['shell']
        df_status = single_job_run(sh, df_status, status_file)

    if is_success(df_status, "index"):
        log_info("Index finished successfully. Proceeding to eqclass.")
    else:
        sh = df_index.iloc[0]['shell'] if not df_index.empty else "unknown"
        log_error(f"Index failed. Check {sh}.stderr")
        return

    # ----------------------
    # Step 2: eqclass
    # ----------------------
    thread_n = int(cfg.getint('xaem', 'eqclass_cpu') / 2)
    count = 1
    while not is_success(df_status, 'eqclass') and count <= 2:
        df_status_eqclass = df_status.query("name == 'eqclass'")
        shell_lst = df_status_eqclass.query("status != 'Success'")['shell'].to_list()

        workQueue = queue.Queue(len(shell_lst))
        for fi in shell_lst:
            workQueue.put(fi)

        threads = []
        actual_threads = min(thread_n, len(shell_lst))
        for i in range(actual_threads):
            thread = MyThread(workQueue, df_status)
            thread.start()
            df_status = thread.df_status
            threads.append(thread)

        for thread in threads:
            thread.join()

        count += 1
        write_status(df_status, status_file)

    if is_success(df_status, "eqclass"):
        success_count = df_status.query("name == 'eqclass' and status == 'Success'").shape[0]
        log_info(f"All eqclass jobs finished successfully ({success_count} jobs). Proceeding to matrix.")
    else:
        error_jobs = df_status.query("name == 'eqclass' and status != 'Success'")
        log_error(f"eqclass failed for {len(error_jobs)} jobs. Check {status_file}.Error for details.")
        return

    # ----------------------
    # Step 3: matrix
    # ----------------------
    df_matrix = df_status[df_status['name'] == 'matrix']
    if not is_success(df_status, "matrix") and not df_matrix.empty:
        sh = df_matrix.iloc[0]['shell']
        df_status = single_job_run(sh, df_status, status_file)

    if is_success(df_status, "matrix"):
        log_info("Matrix step finished. All jobs successfully completed!")
    else:
        sh = df_matrix.iloc[0]['shell'] if not df_matrix.empty else "unknown"
        log_error(f"Matrix step failed. Check {sh}.stderr")

#----------------------------------------------------------------------------
# CLI entry
# ----------------------------------------------------------------------------

@click.command()
@click.option('--verbose', is_flag=True, help='Enable verbose output')
@click.option('-i', '--infile', required=True, type=click.Path(exists=True),
              help='File for sample information (sample name\\tdata source name\\tfq1\\tf2)')
@click.option('--ref', type=click.Choice(['refseq_38', 'gencode_38','pig_110']),
              default='gencode_38', help='Reference transcript')
@click.option('-c', '--config', type=click.Path(exists=True), help='Configuration file')
@click.option('-o', '--outdir', default='./workdir', help='Output directory')
@click.option('--xaem-dir', help='XAEM directory')
@click.option('--xaem-index', help='Pre-built XAEM index')
@click.option('--x-matrix', help='X matrix file')
@click.option('--force', is_flag=True, help='Force to restart all jobs')
def isoquan(verbose, infile, ref, **kwargs):
    log_file = f'{datetime.datetime.now():%Y-%m-%d}.isoquan.info.log'
    if os.path.exists(log_file):
        os.remove(log_file)
    setup_logger(log_file, verbose)
    log_info(f"Project starting. Log file: {log_file}")
    run_isoquan(infile, ref, **kwargs)

