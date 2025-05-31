#!/usr/bin/env python

import os
import logging
import datetime
import configparser
import subprocess
import threading
import queue
from time import asctime
from ..tools import pathfinder, common

#from ..tools.pathfinder import PathFinder
#from ..tools.common import check_file_exists

import pandas as pd
import click

# Configure logging
FORMAT = '%(asctime)s %(message)s'
logger = logging.getLogger(__name__)
binfinder = pathfinder.BinPathFinder('isomap')

class JobStatus(object):
    def __init__(self, sh_nm, df_status):
        self.sh_nm = sh_nm
        self.df_status = df_status
        self.run()

    def run(self):
        logger.warning(f"Running: {self.sh_nm}")
        res = subprocess.run(
            f'bash {self.sh_nm} 1>{self.sh_nm}.stdout 2>{self.sh_nm}.stderr',
            shell=True
        )
        if res.returncode == 0:
            logger.warning(f"Success: {self.sh_nm}")
            self.change_status('Success')
        else:
            logger.warning(f"Error: Please check {self.sh_nm}.stderr")
            self.change_status('Error')
    
    def change_status(self, status_item):
        threadLock.acquire()
        self.df_status.loc[self.df_status['shell'] == self.sh_nm, 'status'] = status_item
        threadLock.release()

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

threadLock = threading.Lock()

def check_fq(fqfile):
    suffixs = ('fq.gz', 'fq', 'fastq', 'fastq.gz', 'fa.gz', 'fa', 'fasta', 'fasta.gz')
    if not fqfile.endswith(suffixs):
        logger.error(f'This is not an fq/fa file, suffix should be {";".join(suffixs)}')
        raise click.BadParameter(f"Invalid file format: {fqfile}")

def read_sampleinfo(infile):
    sample_info = []
    with open(infile) as f:
        for line in f:
            if line.strip() == '' or line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            sample, lib, fq1, fq2 = parts[0:4]
            if not os.path.exists(fq1):
                logger.error(f'{fq1} in {sample} not exists!!')
                raise click.BadParameter(f"File not found: {fq1}")
            elif not os.path.exists(fq2):
                logger.error(f'{fq2} in {sample} not exists!!')
                raise click.BadParameter(f"File not found: {fq2}")
            else:
                sample_info.append([sample, lib, os.path.abspath(fq1), os.path.abspath(fq2)])
    df = pd.DataFrame(sample_info, columns=['sample', 'lib', 'fq1', 'fq2'])
    return df

def index_ref(outdir, config, xaem_dir, refdb, step=1):
    transcript = config.get('xaem', 'transcript_fa') if config.get('xaem', 'transcript_fa') else str(binfinder.find(f'../../data/ref/{refdb}/transcript.fa.gz'))
    common.check_file_exists(
        transcript,
        file_description=f"transcript file is {transcript} for {refdb}",
        logger=logger
    )
    
    cmd = ''
    outfa = f'{outdir}/ref/{os.path.basename(transcript).replace(".gz", "")}'
    index_dir = f'{outdir}/ref/TxIndexer_idx'

    if transcript.endswith('gz'):
        cmd += f'gunzip -c {transcript} > {outfa}\n'
    else:
        cmd += f'ln -fs {transcript} {outfa}\n'
    cmd += f"sed -i 's/|/ /' {outfa}\n"
    cmd += f'{xaem_dir}/bin/TxIndexer -t {outfa} --out {index_dir}\n'
    
    shell_file = f'{outdir}/shell/Step{step}.index_fa.sh'
    os.makedirs(os.path.dirname(shell_file), exist_ok=True)
    with open(shell_file, 'w') as outf:
        outf.write(cmd)
    return index_dir, shell_file

def get_eqclass(df, outdir, xaem_dir, TxIndexer_idx, step=2):
    shell_lst = []
    seqdir = f'{outdir}/seqData'
    os.makedirs(seqdir, exist_ok=True)
    
    for sample, val in df.groupby('sample'):
        cmd = ''
        fq1_lst = list(val['fq1'])
        fq2_lst = list(val['fq2'])
        sample_fq1 = f'{seqdir}/{sample}_1.fq.gz'
        sample_fq2 = f'{seqdir}/{sample}_2.fq.gz'

        if len(fq1_lst) == 1:
            cmd += f'ln -fs {fq1_lst[0]} {sample_fq1}\n'
            cmd += f'ln -fs {fq2_lst[0]} {sample_fq2}\n'
        elif len(fq1_lst) > 1:
            cmd += f'zcat {" ".join(fq1_lst)} | gzip -cf > {sample_fq1} &\n'
            cmd += f'zcat {" ".join(fq2_lst)} | gzip -cf > {sample_fq2} &\n'
            cmd += 'wait\n'
        
        cmd += f"""{xaem_dir}/bin/XAEM \\
        -i {TxIndexer_idx} \\
        -l IU \\
        -1 <(gunzip -c {sample_fq1}) \\
        -2 <(gunzip -c {sample_fq2}) \\
        -p 2 \\
        -o {outdir}/results/{sample}\n"""

        shell_file = f'{outdir}/shell/Step{step}.gen_eqclass_{sample}.sh'
        with open(shell_file, 'w') as outf:
            outf.write(cmd)
        shell_lst.append(shell_file)
    return shell_lst

def count_matrix(outdir, xaem_dir, config, x_matrix, step=3):
    logger.warning(f'Parameter: x_matrix is {x_matrix}')
    resdir = f'{outdir}/results'
    
    
    cmd = f"Rscript {xaem_dir}/R/Create_count_matrix.R workdir={resdir} core=8 design.matrix={x_matrix}\n"
    cmd += f"""Rscript {xaem_dir}/R/AEM_update_X_beta.R \\
        workdir={resdir} \\
        core={config.getint('xaem', 'update_cpu')} \\
        design.matrix={x_matrix} \\
        merge.paralogs={config.getboolean('xaem', 'merge.paralogs')} \\
        isoform.method={config.get('xaem', 'isoform.method')} \\
        remove.ycount={config.getboolean('xaem', 'remove.ycount')}"""
    
    shell_file = f'{outdir}/shell/Step{step}.matrix_samples.sh'
    with open(shell_file, 'w') as outf:
        outf.write(cmd)
    return shell_file

def get_all_shells(outdir, df_sample, config, xaem_dir, refdb, xaem_index=None, x_matrix=None):
    shell_info = []
    step_n = 1
    # Index reference if not provided
    if xaem_index:
        TxIndexer_idx = os.path.abspath(xaem_index)
    else:
        TxIndexer_idx, cmd = index_ref(outdir, config, xaem_dir, refdb, step=step_n)
        shell_info.append([cmd, step_n, 'index'])
        step_n += 1

    # Generate equivalence classes
    eqclass_shells = get_eqclass(df_sample, outdir, xaem_dir, TxIndexer_idx, step=step_n)
    shell_info.extend([[i, step_n, 'eqclass'] for i in eqclass_shells])
    step_n += 1

    # Count matrix and update beta
    if not x_matrix:
        x_matrix = config.get('xaem', 'x_matrix') if config.get('xaem', 'x_matrix') else binfinder.find(f'./resources/ref/{refdb}/X_matrix.RData')
    
    common.check_file_exists(
        x_matrix,
        file_description=f"transcript file is {x_matrix} for {refdb}",
        logger=logger
    )
    
    cmd = count_matrix(outdir, xaem_dir, config, x_matrix, step=step_n)
    shell_info.append([cmd, step_n, 'matrix'])

    # Initialize job status
    df_shell_info = pd.DataFrame(shell_info)
    df_shell_info['status'] = 'Ready'
    df_shell_info.columns = ['shell', 'step', 'name', 'status']
    status_file = f'{outdir}/shell/JOB.Status'
    df_shell_info.to_csv(status_file, index=False, sep='|')
    return status_file

def is_success(df, step_name):
    df_step = df[df['name'] == step_name]
    return df_step.shape[0] == (df_step['status'] == 'Success').sum()

def write_status(df, status_file):
    df.to_csv(status_file, sep='|', index=False)
    df[df['status'] == 'Error'].to_csv(f'{status_file}.Error', sep='|', index=False)

def single_job_run(sh_nm, df, status_file):
    job = JobStatus(sh_nm, df)
    df_status = job.df_status
    write_status(df_status, status_file)
    return df_status

def run_isoquan(infile, ref, config, outdir, xaem_dir, xaem_index, x_matrix, force):
    """Run the full isoform quantification"""
    # Read configuration

    cfg = configparser.ConfigParser()
    if config:
        cfg.read(config, encoding="utf-8")
    else:
        cfg.read(binfinder.find(f'./config.ini'), encoding="utf-8")

    # Set XAEM directory
    if not xaem_dir:
        xaem_dir = cfg.get('xaem', 'xaem_dir') if cfg.get('xaem', 'xaem_dir') else binfinder.find(f'./resources/XAEM/XAEM-binary-0.1.1-cq')
    logger.warning(f'Parameter: xaem_dir is {xaem_dir}')
    if not xaem_dir.exists():
        raise FileNotFoundError(f"XAEM binary not found at {xaem_dir}")
    


    # Create output directories
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(f'{outdir}/seqData', exist_ok=True)
    os.makedirs(f'{outdir}/results', exist_ok=True)
    os.makedirs(f'{outdir}/shell', exist_ok=True)
    os.makedirs(f'{outdir}/ref', exist_ok=True)

    # Read sample information
    df_sample = read_sampleinfo(infile)

    # Get all shell scripts
    status_file = f'{outdir}/shell/JOB.Status'
    if not os.path.exists(status_file) or force:
        logger.warning('Generating all jobs')
        status_file = get_all_shells(
            outdir, df_sample, cfg, xaem_dir, ref, 
            xaem_index=xaem_index, x_matrix=x_matrix
        )
    else:
        logger.warning(f'{status_file} exists. Continuing unfinished jobs. Use --force to restart.')

    # Load job status
    df_status = pd.read_csv(status_file, sep='|')
    status_dict = df_status.groupby('status').size().to_dict()
    logger.warning(f'There are {df_status.shape[0]} jobs, {status_dict}')

    # Step 1: Index reference (if needed)
    if df_status[df_status['name'] == 'index'].shape[0] == 1:
        if not is_success(df_status, 'index'):
            shell_lst = df_status[df_status['name'] == 'index']['shell'].to_list()
            df_status = single_job_run(shell_lst[0], df_status, status_file)

    if is_success(df_status, 'index'):
        logger.warning('Index reference finished, Starting eqclass')
    else:
        logger.warning(f'Index Error, please check {shell_lst[0]}.stderr')

    # Step 2: Process equivalence classes
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

    # Step 3: Generate count matrix
    if is_success(df_status, 'eqclass'):
        df_status_eqclass = df_status.query("name == 'eqclass'")
        success_lst = df_status_eqclass.query("status == 'Success'")['shell'].to_list()
        logger.warning(f'{len(success_lst)} eqclass finished, Starting matrix')
        shell_lst = df_status[df_status['name'] == 'matrix']['shell'].to_list()
        df_status = single_job_run(shell_lst[0], df_status, status_file)
    else:
        error_lst = df_status_eqclass.query("status != 'Error'")['shell'].to_list()
        logger.warning(f'There are {len(error_lst)} errors, please check {status_file}.Error carefully!')

    if is_success(df_status, 'matrix'):
        logger.warning('All jobs finished successfully')
    else:
        logger.warning(f'Matrix not finished. Please check {shell_lst[0]} carefully!')



# 配置日志（保持不变）
FORMAT = '%(asctime)s %(message)s'
logger = logging.getLogger(__name__)

@click.command()
@click.option('--verbose', is_flag=True, help='Enable verbose output')
@click.option('-i', '--infile', required=True, type=click.Path(exists=True),
              help='File for sample information (sample name\\tdata source name\\tfq1\\tfq2)')
@click.option('--ref', type=click.Choice(['refseq_38', 'gencode_38']), 
              default='gencode_38', help='Reference transcript')
@click.option('-c', '--config', type=click.Path(exists=True), 
              help='Configuration file')
@click.option('-o', '--outdir', default='./workdir', 
              help='Output directory')
@click.option('--xaem-dir', help='XAEM directory')
@click.option('--xaem-index', help='Pre-built XAEM index')
@click.option('--x-matrix', help='X matrix file')
@click.option('--force', is_flag=True, help='Force to restart all jobs')
def isoquan(verbose, infile, ref, **kwargs):
    """Isoform quantification"""
    # 设置日志路径
    log_file = f'{datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S")}_info.log'

    # 初始化日志（自定义的 setup_logger 里完成 format 和 level 设置）
    common.setup_logger(log_file, verbose)
    
    logger = logging.getLogger(__name__)
    logger.info(f'Project starting\nLog file: {log_file}')
    
    # 调用核心逻辑
    run_isoquan(infile=infile, ref=ref, **kwargs)

# # 保留直接执行能力
# if __name__ == '__main__':
#     isoquan()