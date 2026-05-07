#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Step S1e: Integration of extra binners for BASALT.

This module provides wrappers around external binning tools (VAMB and
LorBin) and harmonises their outputs into BASALT's internal binset
representation. MetaBinner support was removed in v1.2.x.
"""

from Bio import SeqIO
import subprocess
import sys, os, time
from collections import Counter
from multiprocessing import Pool


def _resolve_vamb_subcmd():
    """Return the VAMB invocation prefix matching the installed major version.

    VAMB 5.x reorganised the CLI into subcommands — default binning is
    ``vamb bin default ...``. VAMB 4.x kept the bare ``vamb ...`` form
    but is pinned to Python 3.11 on bioconda, so within a Python-3.12
    BASALT env users will almost always have v5. Detect at runtime so
    a user-supplied 4.x install in a separate env still works.
    """
    try:
        out = subprocess.check_output(
            ['vamb', '--version'], stderr=subprocess.STDOUT, text=True, timeout=10)
    except (OSError, subprocess.SubprocessError):
        # Fall back to v5 CLI; if vamb is missing entirely the os.system
        # call below will fail loudly with the same error either way.
        return 'vamb bin default'
    out = out.strip()
    # `vamb --version` prints e.g. "5.0.4" or "4.1.3"
    major = out.split('.', 1)[0].lstrip('v').strip()
    return 'vamb' if major == '4' else 'vamb bin default'


def vamb(assembly_file, datasets, num_threads, pwd, QC_software):
    """
    Run VAMB on a single assembly using mapped BAM files.

    Parameters
    ----------
    assembly_file : str
        Assembly FASTA file.
    datasets : dict
        Mapping dataset_id -> [R1, R2] reads used to construct BAMs.
    num_threads : int
        Number of threads to use.
    pwd : str
        Working directory.
    QC_software : {'checkm', 'checkm2'}
        Quality control backend for evaluating resulting bins.
    """
    assembly_num=str(assembly_file).split('_')[0]
    # VAMB requires sorted BAMs (4.x warns on unsorted, 5.x rejects).
    # S1's mapping_paired emits both unsorted (<n>_DNA-<i>.bam) and sorted
    # (<n>_DNA-<i>_sorted.bam) — use the sorted ones.
    for i in range(1, len(datasets)+1):
        if i == 1:
            bam_list=str(assembly_num)+'_DNA-'+str(i)+'_sorted.bam'
        else:
            bam_list+=' '+str(assembly_num)+'_DNA-'+str(i)+'_sorted.bam'

    # VAMB 5.x is required for Python 3.12 (4.x pins to 3.11, conflicts
    # with BASALT's base interpreter). v5 exposes a subcommand layout —
    # `vamb bin default ...` — replacing the legacy bare `vamb ...` CLI.
    vamb_subcmd = _resolve_vamb_subcmd()
    cmd = (vamb_subcmd
           + ' --outdir ' + str(assembly_file) + '_vamb'
           + ' --fasta '  + str(assembly_file)
           + ' --bamfiles '+ str(bam_list)
           + ' --minfasta 500000')
    print('Starting VAMB autobinner: ' + cmd)
    os.system(cmd)
    # os.system('vamb --outdir '+str(assembly_file)+'_100_vamb_genomes --fasta '+str(assembly_file)+' --bamfiles '+str(bam_list)+' -o C')

    vamb_bin_contig, vbn={}, {}
    for line in open(pwd+'/'+str(assembly_file)+'_vamb/clusters.tsv','r'):
        bin_id=str(line).strip().split('\t')[0]
        contig=str(line).strip().split('\t')[1]
        vamb_bin_contig[contig]=bin_id
        vbn[str(assembly_file)+'_100_vamb_genomes.'+str(bin_id)+'.fa']=1

    for record in SeqIO.parse(assembly_file,'fasta'):
        try:
            fmbn=open(str(assembly_file)+'_100_vamb_genomes.'+str(vamb_bin_contig[record.id])+'.fa','a')
            fmbn.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
            fmbn.close()
        except:
            if str(record.id) in vamb_bin_contig.keys():
                fmbn=open(str(assembly_file)+'_100_vamb_genomes.'+str(vamb_bin_contig[record.id])+'.fa','w')
                fmbn.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
                fmbn.close()

    os.system('mkdir '+str(assembly_file)+'_100_vamb_genomes')
    for item in vbn.keys():
        # print(item)
        os.system('mv '+str(item)+' '+pwd+'/'+str(assembly_file)+'_100_vamb_genomes')
    
    os.chdir(pwd+'/'+str(assembly_file)+'_100_vamb_genomes')
    for root, dirs, files in os.walk(pwd+'/'+str(assembly_file)+'_100_vamb_genomes'):
        for file in files:
            bin_len=0
            for record in SeqIO.parse(file,'fasta'):
                bin_len+=len(record.seq)
            if bin_len < 500000:
                os.system('rm '+str(file))
    os.chdir(pwd)

    if QC_software == 'checkm':
        os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(assembly_file)+'_100_vamb_genomes '+str(assembly_file)+'_100_vamb_checkm')
    elif QC_software == 'checkm2':
        os.system('checkm2 predict -t '+str(num_threads)+' -i '+str(assembly_file)+'_100_vamb_genomes  -x fa -o '+str(assembly_file)+'_100_vamb_checkm')
    # os.system('rm *.seed *.out *.err *.nto *.gff *.ffn *.faa *.ndb *.njs *.not *.ntf')
    os.system('rm -rf '+str(assembly_file)+'_vamb')


def lorbin(assembly_file, datasets, num_threads, pwd, QC_software):
    """
    Run LorBin on a single assembly using mapped BAM files.

    Parameters
    ----------
    assembly_file : str
        Assembly FASTA file.
    datasets : dict
        Mapping dataset_id -> [R1, R2] reads used to construct BAMs.
    num_threads : int
        Number of threads to use.
    pwd : str
        Working directory.
    QC_software : {'checkm', 'checkm2'}
        Quality control backend for evaluating resulting bins.
    """
    assembly_num=str(assembly_file).split('_')[0]
    # LorBin requires sorted+indexed BAMs (it calls samtools internally).
    # Use the *_sorted.bam set produced by S1's mapping_paired.
    for i in range(1, len(datasets)+1):
        if i == 1:
            bam_list=str(assembly_num)+'_DNA-'+str(i)+'_sorted.bam'
        else:
            bam_list+=' '+str(assembly_num)+'_DNA-'+str(i)+'_sorted.bam'

    # Lowercase folder name to match the convention used by every other
    # binner output (`_100_metabat_genomes`, `_100_vamb_genomes`, ...) and
    # to match the substring matching in s2/s6/s7 filename detection.
    lorbin_genome = str(assembly_file) + '_100_lorbin_genomes'

    bam_list = str(bam_list).split()
    bam_abs_list = [os.path.join(pwd, b) for b in bam_list]
    bam_abs = ' '.join(bam_abs_list)

    # 多样本时加 --multi
    multi_flag = ''
    if len(bam_list) > 1:
        multi_flag = ' --multi'

    cmd = (
        'LorBin bin'
        ' -o ' + lorbin_genome +
        ' -fa ' + str(assembly_file) +
        ' -b ' + bam_abs +
        ' --num_process ' + str(num_threads) +
        multi_flag
    )

    print('Starting LorBin autobinner')
    print(cmd)
    os.system(cmd)
    print('LorBin binning finished for ' + lorbin_genome)

    # Normalise bin filenames to: <assembly_file>_100_lorbin_genomes.N.fa
    bin_index = 0
    outdir = os.path.join(pwd, lorbin_genome)
    for root, dirs, files in os.walk(outdir):
        for file in files:
            if file.endswith('.fa') or file.endswith('.fna') or file.endswith('.fasta'):
                bin_index += 1
                old_path = os.path.join(root, file)
                new_name = lorbin_genome + '.' + str(bin_index) + '.fa'
                new_path = os.path.join(outdir, new_name)
                if old_path != new_path:
                    os.system('mv ' + old_path + ' ' + new_path)
    os.chdir(pwd)

    if QC_software == 'checkm':
        os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+lorbin_genome+' '+str(assembly_file)+'_100_lorbin_checkm')
    elif QC_software == 'checkm2':
        os.system('checkm2 predict -t '+str(num_threads)+' -i '+lorbin_genome+'  -x fa -o '+str(assembly_file)+'_100_lorbin_checkm')



def extra_binner(binner, datasets, assembly_file, depth_file, num_threads, ram, pwd, QC_software):
    """
    Dispatcher function to select and run the appropriate extra binner.

    Supports VAMB ('v') and LorBin ('l'), and returns a list of new
    binset folder names produced by the selected tool.
    """
    extra_bin_folder =[]

    # VAMB
    if binner == 'v':
        vamb(assembly_file, datasets, num_threads, pwd, QC_software)
        extra_bin_folder.append(str(assembly_file) + '_100_vamb_genomes')

    elif binner == 'l':
        lorbin(assembly_file, datasets, num_threads, pwd, QC_software)
        extra_bin_folder.append(str(assembly_file) + '_100_lorbin_genomes')

    # Cleanup of common temporary files generated by these tools
    os.system('rm *.seed *.out *.err *.nto *.gff *.ffn *.faa *.ndb *.njs *.not *.ntf 2>/dev/null')

    return extra_bin_folder

if __name__ == '__main__':
    num_threads=20
    ram=250
    pwd=os.getcwd()
    assembly_file='1_assembly_sample1.fa'
    depth_file='1_assembly.depth.txt'
    datasets={'1':['sample1.R1.fq','sample1.R2.fq'], '2':['sample2.R1.fq','sample2.R2.fq']}
    binner='v' ### 'v': vamb; 'l': lorbin
    QC_software='checkm2' ### checkm or checkm2
    extra_binner(binner, datasets, assembly_file, depth_file, num_threads, ram, pwd, QC_software)
