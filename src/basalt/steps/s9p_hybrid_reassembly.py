#!/usr/bin/env python

"""
Step S9p: Hybrid reassembly using both short and long reads.

Unified entry for CheckM2 and CheckM backends. Call ``set_qc_backend(name)``
before invoking ``hybrid_re_assembly_main`` to select the QC tool.
"""

from Bio import SeqIO
import sys, os, threading, copy, math
from multiprocessing import Pool

from basalt.qc_backend import get_backend, strip_fasta_extension


# Module-level QC backend. Set by the caller (BASALT_main_autobinning / BASALT_main_refinement / BASALT_main_re_assembly)
# via ``set_qc_backend`` before running the pipeline so that all helpers in
# this module know whether to invoke CheckM2 or CheckM.
_BACKEND = None


def set_qc_backend(qc_software):
    """Configure the QC backend used by all helpers in this module."""
    global _BACKEND
    _BACKEND = get_backend(qc_software)


def _require_backend():
    if _BACKEND is None:
        raise RuntimeError(
            'S9p_Hybrid_Reassembly: QC backend not set. '
            "Call set_qc_backend('checkm2' or 'checkm') first."
        )
    return _BACKEND


def hybrid_parse_checkm(checkm_containing_folder, pwd):
    """
    Parse QC reports produced by either CheckM2 or CheckM for S9p.

    Always call with the top-level output folder. The backend walks the
    tree, so CheckM's ``storage/bin_stats_ext.tsv`` is found without the
    caller having to append ``/storage``.
    """
    backend = _require_backend()
    raw = backend.parse_results(os.path.join(pwd, checkm_containing_folder))
    return {strip_fasta_extension(bin_id): metrics for bin_id, metrics in raw.items()}


def assembly_mul(bins_seq_folder, bin_seq, item, reassembly_bin_folder,
                 pwd, num_threads, ram, sensitive='sensitive'):
    """
    Per-bin short-read reassembly (CheckM2 branch). Always runs MEGAHIT
    (--presets meta-sensitive); also runs SPAdes when sensitive=='more-sensitive'.
    """
    try:
        fx=open('Basalt_log.txt','a')
    except:
        fx=open('Basalt_log.txt','w')

    megahit_out=str(item)+'_megahit_reassembly'
    fx.write('Assembling '+str(item)+' using MEGAHIT (--presets meta-sensitive)'+'\n')
    os.system('rm -rf '+megahit_out)
    os.system('megahit --presets meta-sensitive -1 '+pwd+'/'+bins_seq_folder+'/'+str(item)+'_seq_R1.fq -2 '+pwd+'/'+bins_seq_folder+'/'+str(item)+'_seq_R2.fq -o '+megahit_out+' --min-contig-len 1000 -t '+str(num_threads))
    mh_n=0
    if os.path.isfile(megahit_out+'/final.contigs.fa'):
        f=open(str(item)+'_megahit_re-assembly_contigs.fa','w')
        for record in SeqIO.parse(megahit_out+'/final.contigs.fa', 'fasta'):
            if len(record.seq) >= 1000:
                mh_n+=1
                f.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
        f.close()
        if mh_n >= 2:
            os.system('mv '+str(item)+'_megahit_re-assembly_contigs.fa '+str(reassembly_bin_folder)+'/'+str(item)+'_megahit_re-assembly_contigs.fa')
    os.system('rm -rf '+megahit_out)

    if sensitive == 'more-sensitive':
        fx.write('Assembling '+str(item)+' using SPAdes'+'\n')
        os.system('spades.py -1 '+pwd+'/'+bins_seq_folder+'/'+str(item)+'_seq_R1.fq -2 '+pwd+'/'+bins_seq_folder+'/'+str(item)+'_seq_R2.fq -o '+str(item)+'_spades_reassembly --careful -t '+str(num_threads)+' -m '+str(ram))
        os.system('mv '+pwd+'/'+str(item)+'_spades_reassembly/contigs.fasta '+str(item)+'_contigs.fasta')

        xxxx=0
        if os.path.isfile(str(item)+'_contigs.fasta') == True:
            f=open(str(item)+'_SPAdes_re-assembly_contigs.fa','w')
            for record in SeqIO.parse(str(item)+'_contigs.fasta', 'fasta'):
                if len(record.seq) >= 1000:
                    xxxx+=1
                    f.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
            f.close()

            os.system('mv '+pwd+'/'+str(item)+'_spades_reassembly/corrected/'+str(item)+'_seq_R1* '+pwd+'/SPAdes_corrected_reads/'+str(item)+'_seq_R1.fq.gz')
            os.system('mv '+pwd+'/'+str(item)+'_spades_reassembly/corrected/'+str(item)+'_seq_R2* '+pwd+'/SPAdes_corrected_reads/'+str(item)+'_seq_R2.fq.gz')
        else:
            os.system('spades.py -1 '+pwd+'/'+bins_seq_folder+'/'+str(item)+'_seq_R1.fq -2 '+pwd+'/'+bins_seq_folder+'/'+str(item)+'_seq_R2.fq -o '+str(item)+'_spades_reassembly -t '+str(num_threads)+' -m '+str(ram))
            os.system('mv '+pwd+'/'+str(item)+'_spades_reassembly/contigs.fasta '+str(item)+'_contigs.fasta')
            if os.path.isfile(str(item)+'_contigs.fasta') == True:
                f=open(str(item)+'_SPAdes_re-assembly_contigs.fa','w')
                for record in SeqIO.parse(str(item)+'_contigs.fasta', 'fasta'):
                    if len(record.seq) >= 1000:
                        xxxx+=1
                        f.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
                f.close()
            os.system('mv '+pwd+'/'+str(item)+'_spades_reassembly/corrected/'+str(item)+'_seq_R1* '+pwd+'/SPAdes_corrected_reads/'+str(item)+'_seq_R1.fq.gz')
            os.system('mv '+pwd+'/'+str(item)+'_spades_reassembly/corrected/'+str(item)+'_seq_R2* '+pwd+'/SPAdes_corrected_reads/'+str(item)+'_seq_R2.fq.gz')

        os.system('rm '+str(item)+'_contigs.fasta')
        if xxxx >= 2:
            os.system('mv '+str(item)+'_SPAdes_re-assembly_contigs.fa '+str(reassembly_bin_folder)+'/'+str(item)+'_SPAdes_re-assembly_contigs.fa')
        os.system('rm -rf '+str(item)+'_spades_reassembly')

    fx.close()

def SR_reassembly(bin_seq, reassembly_bin_folder, num_threads,
                  bins_seq_folder, long_read, ram, pwd, sensitive='sensitive'):
    """
    Run short-read-only hybrid reassembly workflow (CheckM2 branch).

    Parameters
    ----------
    bin_seq : dict
        Mapping ``bin_id -> [r1_fastq, r2_fastq]`` with bin-specific reads.
    reassembly_bin_folder : str
        Folder where reassembled bin FASTA files will be written.
    num_threads : int
        Number of threads for SPAdes.
    bins_seq_folder : str
        Folder containing bin-specific short-read FASTQ files.
    long_read : list of str
        Long-read datasets (used indirectly when scheduling hybrid runs).
    ram : int
        Maximum RAM (in GB) available to assemblers.
    pwd : str
        Working directory path.

    Returns
    -------
    dict
        Mapping ``bin_id -> dict(checkm_metrics)`` for reassembled bins.
    """
    num_project=1
    if num_threads >= 40:
        if num_threads < 60:
            num_project=2
        else:
            num_project=math.ceil(num_threads/30)
    if ram >= 64:
        num_project2=math.ceil(ram/55)
        if num_project2 < num_project:
            num_project=num_project2
    
    os.system('mkdir SPAdes_corrected_reads')
    # for item in bin_seq.keys():
    #     os.system('spades.py -1 '+str(pwd)+'/'+str(bins_seq_folder)+'/'+str(bin_seq[item][0])+' -2 '+str(pwd)+'/'+str(bins_seq_folder)+'/'+str(bin_seq[item][1])+' -o '+str(item)+'_SPAdes_corrected_reads --only-error-correction -t '+str(num_threads)+' -m '+str(ram))
    #     os.system('mv '+pwd+'/'+str(item)+'_SPAdes_corrected_reads/corrected/'+str(item)+'_seq_R1* '+pwd+'/SPAdes_corrected_reads/'+str(item)+'_seq_R1.fq.gz')
    #     os.system('mv '+pwd+'/'+str(item)+'_SPAdes_corrected_reads/corrected/'+str(item)+'_seq_R2* '+pwd+'/SPAdes_corrected_reads/'+str(item)+'_seq_R2.fq.gz')
    #     os.system('rm -rf '+str(item)+'_SPAdes_corrected_reads')

    t_p_p=math.ceil(num_threads/num_project)
    print('Reassembly parallelism: '+str(num_project)+' bins assembled in parallel x '+str(t_p_p)+' threads/bin '
          '(--threads='+str(num_threads)+', --ram='+str(ram)+'G; assemblers will report '+str(t_p_p)+' threads each by design)')
    pool=Pool(processes=num_project)
    for item in bin_seq.keys():
        print('Reassembling '+str(item))
        pool.apply_async(assembly_mul, args=(bins_seq_folder, bin_seq, item, reassembly_bin_folder, pwd, t_p_p, ram, sensitive))
        # assembly_mul('SPAdes_corrected_reads', bins_seq_folder, bin_seq, item, reassembly_bin_folder, pwd, t_p_p, ram)
    pool.close()
    pool.join()

def hybrid_bin_comparison(paired_bins, bin_checkm):
    # pwd=os.getcwd()
    f=open('Reassembled_bins_comparison.txt','w')
    best_bin, best_bin_checkm={}, {}
    for item in paired_bins.keys():
        best_bin_checkm_name=None
        for item2 in paired_bins[item]:
            if '_polished' in item2:
                best_bin_checkm_name_list=item2.split('.')
                best_bin_checkm_name_list.remove(best_bin_checkm_name_list[-1])
                best_bin_checkm_name='.'.join(best_bin_checkm_name_list)
                # best_bin_cpn=bin_checkm[best_bin_checkm_name]['Completeness']
                # best_bin_ctn=bin_checkm[best_bin_checkm_name]['Contamination']
                # best_bin_ml=bin_checkm[best_bin_checkm_name]['Mean scaffold length']
                f.write(str(item2)+'\t'+str(bin_checkm[best_bin_checkm_name])+'\n')

        for item2 in paired_bins[item]:
            if '_polished' not in item2:
                reass_bin_checkm_name_list=item2.split('.')
                reass_bin_checkm_name_list.remove(reass_bin_checkm_name_list[-1])
                reass_bin_checkm_name='.'.join(reass_bin_checkm_name_list)
                f.write(str(item2)+'\t'+str(bin_checkm[reass_bin_checkm_name])+'\n')
                if best_bin_checkm_name is None:
                    # No polished baseline for this bin; nothing to compare
                    # against, so adopt the reassembly entry as best.
                    best_bin_checkm_name=reass_bin_checkm_name
                    continue
                best_bin_cpn=bin_checkm[best_bin_checkm_name]['Completeness']
                best_bin_ctn=bin_checkm[best_bin_checkm_name]['Contamination']
                best_bin_ml=bin_checkm[best_bin_checkm_name]['contig_size']
                reass_bin_cpn=bin_checkm[reass_bin_checkm_name]['Completeness']
                reass_bin_ctn=bin_checkm[reass_bin_checkm_name]['Contamination']
                reass_bin_ml=bin_checkm[reass_bin_checkm_name]['contig_size']
        
                delta_cpn_ctn_bestbin=float(best_bin_cpn)-5*float(best_bin_ctn)
                delta_cpn_ctn_reass_bin=float(reass_bin_cpn)-float(5*reass_bin_ctn)

                delta_cpn_ctn_bestbin_x=float(best_bin_cpn)-float(best_bin_ctn)
                delta_cpn_ctn_reass_bin_x=float(reass_bin_cpn)-float(reass_bin_ctn)

                if '_hybird_' in reass_bin_checkm_name:
                    if delta_cpn_ctn_bestbin > delta_cpn_ctn_reass_bin:
                        if delta_cpn_ctn_bestbin >= 40:
                            best_bin_checkm_name=best_bin_checkm_name
                        else:
                            if delta_cpn_ctn_bestbin_x >= delta_cpn_ctn_reass_bin_x:
                                best_bin_checkm_name=best_bin_checkm_name
                            elif delta_cpn_ctn_bestbin_x < delta_cpn_ctn_reass_bin_x: 
                                best_bin_checkm_name=reass_bin_checkm_name
                    elif delta_cpn_ctn_bestbin < delta_cpn_ctn_reass_bin:
                        if delta_cpn_ctn_reass_bin >= 40:
                            best_bin_checkm_name=reass_bin_checkm_name
                        else:
                            if delta_cpn_ctn_bestbin_x >= delta_cpn_ctn_reass_bin_x:
                                best_bin_checkm_name=best_bin_checkm_name
                            elif delta_cpn_ctn_bestbin_x < delta_cpn_ctn_reass_bin_x: 
                                best_bin_checkm_name=reass_bin_checkm_name
                    elif delta_cpn_ctn_bestbin == delta_cpn_ctn_reass_bin:
                        if reass_bin_ml > best_bin_ml:
                            best_bin_checkm_name=reass_bin_checkm_name
                        else:
                            best_bin_checkm_name=best_bin_checkm_name
                    else:
                        continue
                else:
                    delta_idba = delta_cpn_ctn_reass_bin-delta_cpn_ctn_bestbin
                    if delta_idba < 3:
                        if delta_cpn_ctn_bestbin >= 40:
                            best_bin_checkm_name=best_bin_checkm_name
                        else:
                            if delta_cpn_ctn_bestbin_x >= delta_cpn_ctn_reass_bin_x:
                                best_bin_checkm_name=best_bin_checkm_name
                            elif delta_cpn_ctn_bestbin_x < delta_cpn_ctn_reass_bin_x: 
                                best_bin_checkm_name=reass_bin_checkm_name
                    elif delta_idba >= 3 and delta_idba < 6:
                        ratio=reass_bin_ml/best_bin_ml
                        if ratio >= 0.5:
                            if delta_cpn_ctn_reass_bin >= 40:
                                best_bin_checkm_name=reass_bin_checkm_name
                            else:
                                if delta_cpn_ctn_bestbin_x >= delta_cpn_ctn_reass_bin_x:
                                    best_bin_checkm_name=best_bin_checkm_name
                                elif delta_cpn_ctn_bestbin_x < delta_cpn_ctn_reass_bin_x: 
                                    best_bin_checkm_name=reass_bin_checkm_name
                        else:
                            best_bin_checkm_name=best_bin_checkm_name
                    elif delta_idba >= 6:
                        if delta_cpn_ctn_reass_bin >= 40:
                            best_bin_checkm_name=reass_bin_checkm_name
                        else:
                            if delta_cpn_ctn_bestbin_x >= delta_cpn_ctn_reass_bin_x:
                                best_bin_checkm_name=best_bin_checkm_name
                            elif delta_cpn_ctn_bestbin_x < delta_cpn_ctn_reass_bin_x: 
                                best_bin_checkm_name=reass_bin_checkm_name 
                    else:
                        continue

        if best_bin_checkm_name is None:
            print('Skipping '+str(item)+': no polished or reassembly entry found')
            continue
        best_bin[best_bin_checkm_name+'.fa']=best_bin_checkm_name
        best_bin_checkm[best_bin_checkm_name]=bin_checkm[best_bin_checkm_name].copy()
    f.close()
    return best_bin, best_bin_checkm

def hybrid_assembly_mul(sr_folder, bin_seq, item, bin_lr_reads, lr_folder,
                        reassembly_bin_folder, pwd, num_threads, ram):
    """
    Assemble hybrid (short + long read) data for a single bin using SPAdes
    (CheckM2 branch).

    Parameters
    ----------
    sr_folder : str
        Folder containing short-read FASTQ files used for hybrid assembly.
    bin_seq : dict
        Mapping ``bin_id -> [r1_fastq, r2_fastq]`` with bin-specific reads.
    item : str
        Bin identifier to assemble.
    bin_lr_reads : dict
        Mapping ``bin_id -> long_read_fastq`` with bin-specific long reads.
    lr_folder : str
        Folder containing long-read FASTQ files.
    reassembly_bin_folder : str
        Folder where hybrid reassembled bin FASTA files will be written.
    pwd : str
        Working directory path.
    num_threads : int
        Number of threads to use for SPAdes.
    ram : int
        Maximum RAM (in GB) available to SPAdes.

    Returns
    -------
    None
        Writes hybrid reassembled contig FASTA files into
        ``reassembly_bin_folder`` and uses ``SPAdes_corrected_reads``
        to store corrected short reads.
    """
    try:
        os.system('gzip -d '+pwd+'/SPAdes_corrected_reads/'+str(item)+'_seq_R1.fq.gz')
        os.system('gzip -d '+pwd+'/SPAdes_corrected_reads/'+str(item)+'_seq_R2.fq.gz')
    except:
        print(str(item)+' corrected reads does not exist')

    # try:
    #     os.system('mv '+pwd+'/'+bins_seq_folder+'/'+str(item)+'_seq_R1.fq.gz '+pwd+'/'+bins_seq_folder+'/'+str(item)+'_seq_R2.fq.gz '+pwd)
    #     os.system('gzip -d '+str(item)+'_seq_R1.fq.gz')
    #     os.system('gzip -d '+str(item)+'_seq_R2.fq.gz')
    # except:
    #     print(str(item)+' corrected reads does not exist')


    if os.path.isfile(pwd+'/'+sr_folder+'/'+str(item)+'_seq_R1.fq') == True:
        os.system('spades.py -1 '+pwd+'/SPAdes_corrected_reads/'+str(item)+'_seq_R1.fq -2 '+pwd+'/SPAdes_corrected_reads/'+str(item)+'_seq_R2.fq --nanopore '+pwd+'/'+lr_folder+'/'+str(bin_lr_reads[item])+' -o '+str(item)+'_spades_hybrid_reassembly --careful --only-assembler -t '+str(num_threads)+' -m '+str(ram))
        # os.system('rm '+str(item)+'_seq_R1.fq '+str(item)+'_seq_R2.fq')
    else:
        os.system('spades.py -1 '+pwd+'/'+sr_folder+'/'+str(bin_seq[item][0])+' -2 '+pwd+'/'+sr_folder+'/'+str(bin_seq[item][1])+' --nanopore '+pwd+'/'+lr_folder+'/'+str(bin_lr_reads[item])+' -o '+str(item)+'_spades_hybrid_reassembly --careful -t '+str(num_threads)+' -m '+str(ram))
    os.system('mv '+str(pwd)+'/'+str(item)+'_spades_hybrid_reassembly/contigs.fasta '+str(item)+'_contigs.fasta')
    # os.chdir(str(pwd)+'/'+str(item)+'_spades_hybrid_reassembly')
    xxxx=0
    if os.path.isfile(str(item)+'_contigs.fasta') == True:
        f=open(str(item)+'_SPAdes_hybrid_re-assembly_contigs.fa','w')
        for record in SeqIO.parse(str(item)+'_contigs.fasta', 'fasta'):
            if len(record.seq) >= 1000:
                xxxx+=1
                f.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
        f.close()
        os.system('rm '+str(item)+'_contigs.fasta')
    # os.chdir(str(pwd))
    
    if xxxx >= 1:
        os.system('mv '+str(item)+'_SPAdes_hybrid_re-assembly_contigs.fa '+str(reassembly_bin_folder)+'/'+str(item)+'_SPAdes_hybrid_re-assembly_contigs.fa')
        # os.system('mv '+str(item)+'_spades_hybrid_reassembly/'+str(item)+'_SPAdes_hybrid_re-assembly_contigs.fa '+str(reassembly_bin_folder)+'/'+str(item)+'_SPAdes_hybrid_re-assembly_contigs.fa')
    os.system('rm -rf '+str(item)+'_spades_hybrid_reassembly')

def hybrid_re_assembly_main(binset_folder, sr_folder, lr_folder,
                            ram, num_threads, sensitive='sensitive'):
    """
    Entry point for S9p hybrid reassembly (CheckM2 branch).

    Parameters
    ----------
    binset_folder : str
        Folder containing the starting binset.
    sr_folder : str
        Folder holding short-read FASTQ files used in hybrid assembly.
    lr_folder : str
        Folder holding long-read FASTQ files.
    ram : int
        Maximum RAM (in GB) available to assemblers.
    num_threads : int
        Number of threads to use for mapping and assembly.

    Returns
    -------
    None
        Coordinates short-read assembly, hybrid assembly, CheckM-based
        bin comparison, and writes results into the hybrid reassembly
        folders on disk.
    """
    pwd=os.getcwd()
    bin_checkm={}
    try:
        f=open('Hybrid_re-assembly_status.txt','a')
    except:
        f=open('Hybrid_re-assembly_status.txt','w')
    f.close()

    reassembly_bin_folder=binset_folder+'_re-assembly_binset'
    try:
        os.mkdir(reassembly_bin_folder)
    except:
        print(reassembly_bin_folder+' exists')

    bin_seq={}
    os.chdir(pwd+'/'+binset_folder)
    for root, dirs, files in os.walk(pwd+'/'+binset_folder):
        for file in files:
            if '_mag_polished.fa' in file:
                bin_id=str(file).split('_mag_polished.fa')[0]
                bin_seq[bin_id]=[]
                os.system('cp '+str(file)+' '+pwd+'/'+reassembly_bin_folder)

    os.chdir(pwd+'/'+sr_folder)
    for root, dirs, files in os.walk(pwd+'/'+sr_folder):
        for file in files:
            if '_seq_R1.fq' in file:
                bin_id=str(file).split('_seq_R1.fq')[0]
                bin_seq[bin_id].append(file)
            elif '_seq_R2.fq' in file:
                bin_id=str(file).split('_seq_R2.fq')[0]
                bin_seq[bin_id].append(file)
    os.chdir(pwd)
    
    x=0
    for line in open('Hybrid_re-assembly_status.txt','r'):
        if 'Short-read assembly done!' in line:
            x=1

    if x == 0:
        SR_reassembly(bin_seq, reassembly_bin_folder, num_threads, sr_folder, lr_folder, ram, pwd, sensitive=sensitive)
        f=open('Hybrid_re-assembly_status.txt','a')
        f.write('Short-read assembly done!'+'\n')
        f.close()

    if lr_folder != '':
        x=0
        for line in open('Hybrid_re-assembly_status.txt','r'):
            if 'Hybrid-assembly done!' in line:
                x=1

        if x == 0:
            bin_lr_reads={}
            os.chdir(pwd+'/'+lr_folder)
            for root, dirs, files in os.walk(pwd+'/'+lr_folder):
                for file in files:
                    if '_lr.fq' in file:
                        bin_id=str(file).split('_lr.fq')[0]
                        bin_lr_reads[bin_id]=file
            os.chdir(pwd)

            ####
            num_project=1
            if num_threads >= 40:
                if num_threads < 60:
                    num_project=2
                else:
                    num_project=math.ceil(num_threads/30)
            if ram >= 64:
                num_project2=math.ceil(ram/55)
                if num_project2 < num_project:
                    num_project=num_project2
            
            t_p_p=math.ceil(num_threads/num_project)
            print('Hybrid reassembly parallelism: '+str(num_project)+' bins assembled in parallel x '+str(t_p_p)+' threads/bin '
                  '(--threads='+str(num_threads)+', --ram='+str(ram)+'G; assemblers will report '+str(t_p_p)+' threads each by design)')
            pool=Pool(processes=num_project)
            for item in bin_lr_reads.keys():
                print('Reassembling '+str(item))
                pool.apply_async(hybrid_assembly_mul, args=(sr_folder, bin_seq, item, bin_lr_reads, lr_folder, reassembly_bin_folder, pwd, t_p_p, ram))
            pool.close()
            pool.join()

            f=open('Hybrid_re-assembly_status.txt','a')
            f.write('Hybrid-assembly done!'+'\n')
            f.close()

    x=0
    for line in open('Hybrid_re-assembly_status.txt','r'):
        if 'Re-assembly quality evaluation done!' in line:
            x=1
    
    if x == 0:
        backend = _require_backend()
        os.system(backend.build_cmd(
            num_threads,
            str(binset_folder) + '_re-assembly_binset',
            str(binset_folder) + '_re-assembly_binset_checkm',
            ext='fa',
        ))
        f=open('Hybrid_re-assembly_status.txt','a')
        f.write('Re-assembly quality evaluation done!'+'\n')
        f.close()

    reassembly_bin_checkm=hybrid_parse_checkm(binset_folder+'_re-assembly_binset_checkm',pwd)
    bin_checkm.update(reassembly_bin_checkm)

    x=0
    for line in open('Hybrid_re-assembly_status.txt','r'):
        if 'Bin selecting done!' in line:
            x=1
    
    if x == 0:
        paired_bins, best_bin, best_bin_checkm={}, {}, {}
        os.chdir(pwd+'/'+binset_folder+'_re-assembly_binset')
        for root, dirs, files in os.walk(pwd+'/'+binset_folder+'_re-assembly_binset'):
            for file in files:
                if '_SPAdes_re-assembly_contigs.fa' in file:
                    n_re=0
                    for line in open(file,'r'):
                        n_re+=1
                        if n_re == 2:
                            break
                    if n_re == 2:
                        o_bin=str(file).split('_SPAdes_re-assembly_contigs.fa')[0]+'.fa'
                        if o_bin not in paired_bins.keys():
                            paired_bins[o_bin]=[]
                            paired_bins[o_bin].append(str(file))
                        else:
                            paired_bins[o_bin].append(str(file))
                elif '_SPAdes_hybrid_re-assembly_contigs.fa' in file:
                    n_re=0
                    for line in open(file,'r'):
                        n_re+=1
                        if n_re == 2:
                            break
                    if n_re == 2:
                        o_bin=str(file).split('_SPAdes_hybrid_re-assembly_contigs.fa')[0]+'.fa'
                        if o_bin not in paired_bins.keys():
                            paired_bins[o_bin]=[]
                            paired_bins[o_bin].append(str(file))
                        else:
                            paired_bins[o_bin].append(str(file))
                elif '_IDBA_re-assembly_contigs.fa' in file:
                    n_re=0
                    for line in open(file,'r'):
                        n_re+=1
                        if n_re == 2:
                            break
                    if n_re == 2:
                        o_bin=str(file).split('_IDBA_re-assembly_contigs.fa')[0]+'.fa'
                        if o_bin not in paired_bins.keys():
                            paired_bins[o_bin]=[]
                            paired_bins[o_bin].append(str(file))
                        else:
                            paired_bins[o_bin].append(str(file))
                elif '_mag_polished.fa' in file:
                    n_re=0
                    for line in open(file,'r'):
                        n_re+=1
                        if n_re == 2:
                            break
                    if n_re == 2:
                        o_bin=str(file).split('_mag_polished.fa')[0]+'.fa'
                        if o_bin not in paired_bins.keys():
                            paired_bins[o_bin]=[]
                            paired_bins[o_bin].append(str(file))
                        else:
                            paired_bins[o_bin].append(str(file))
                else:
                    o_bin=str(file).split('.')[0]
                    if o_bin not in paired_bins.keys():
                        paired_bins[o_bin]=[]
                        paired_bins[o_bin].append(str(file))
                    else:
                        paired_bins[o_bin].append(str(file))
        os.chdir(pwd)

        best_bin_after_c=hybrid_bin_comparison(paired_bins, bin_checkm)
        best_bin=best_bin_after_c[0].copy()
        best_bin_checkm=best_bin_after_c[1].copy()

        try:
            os.mkdir(binset_folder+'_re-assembly')
        except:
            print(binset_folder+'_re-assembly exists')
        
        selected_bins={}
        for item in best_bin.keys():
            os.system('cp '+pwd+'/'+reassembly_bin_folder+'/'+item+' '+pwd+'/'+binset_folder+'_re-assembly')

        os.chdir(pwd+'/'+binset_folder+'_re-assembly')
        backend = _require_backend()
        f=open('Best_binset_after_re-assembly_quality_report.tsv','w')
        f.write('Bin_ID\tGenome_size\tCompleteness\tContamination\t' + backend.contig_size_key + '\n')

        for item in best_bin_checkm.keys():
            f.write(
                str(item) + '\t'
                + str(best_bin_checkm[item]['Genome size']) + '\t'
                + str(best_bin_checkm[item]['Completeness']) + '\t'
                + str(best_bin_checkm[item]['Contamination']) + '\t'
                + str(best_bin_checkm[item]['contig_size']) + '\n'
            )
        f.close()
        os.chdir(pwd)

        f=open('Hybrid_re-assembly_status.txt','a')
        f.write('Bin selecting done!'+'\n')
        f.close()
    print('Re-assembly done!')

if __name__ == '__main__':
    binset_folder='BestBinset_outlier_refined_filtrated_retrieved_MAGs_polished'
    sr_folder='BestBinset_outlier_refined_filtrated_retrieved_polished_sr_bins_seq'
    lr_folder='BestBinset_outlier_refined_filtrated_retrieved_long_read'
    num_threads=60
    ram=250
    set_qc_backend('checkm2')
    hybrid_re_assembly_main(binset_folder, sr_folder, lr_folder, ram, num_threads)
