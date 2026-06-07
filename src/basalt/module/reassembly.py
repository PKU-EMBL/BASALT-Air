#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Re-assembly entry point for the BASALT-Air pipeline.

This module orchestrates OLC-based reassembly, hybrid reassembly, and
final polishing of bins when long reads or HiFi datasets are provided.

Supports both CheckM2 and CheckM via the ``QC_software`` parameter.
"""

import time
import sys
import os
from Bio import SeqIO
from basalt.steps.s4_multiple_assembly_comparator import *
from basalt.steps.s4_multiple_assembly_comparator import set_qc_backend as _set_s4_backend
from basalt.steps.s8_olc import *
from basalt.steps.s8_olc import set_qc_backend as _set_s8_backend
from basalt.steps.s9_reassembly import *
from basalt.steps.s9_reassembly import set_qc_backend as _set_s9_backend
from basalt.steps.s9p_hybrid_reassembly import *
from basalt.steps.s9p_hybrid_reassembly import set_qc_backend as _set_s9p_backend
from basalt.steps.s10_olc import *
from basalt.steps.s10_olc import set_qc_backend as _set_s10_backend
from basalt.core.cleanup import *
from glob import glob  # keep last: cleanup re-exports `glob` (module) via wildcard


def _folder_has_fasta(path):
    if not os.path.isdir(path):
        return False
    for _, _, files in os.walk(path):
        for f in files:
            ext = f.rsplit('.', 1)[-1].lower()
            if ext in ('fa', 'fna', 'fasta'):
                return True
    return False


def _resolve_pre_olc_bin_folder(pwd, lr_list):
    """Pick the bin folder fed into the OLC step when the checkpoint is empty.

    Mirrors the historical inline logic: prefer the polished MAGs when long
    reads are present, otherwise fall back through the retrieved-retrieved
    and retrieved folders.
    """
    candidates = []
    if len(lr_list) != 0:
        candidates.append('BestBinset_outlier_refined_filtrated_retrieved_MAGs_polished')
    candidates.append('BestBinset_outlier_refined_filtrated_retrieved_retrieved')
    candidates.append('BestBinset_outlier_refined_filtrated_retrieved')
    for c in candidates:
        if _folder_has_fasta(os.path.join(pwd, c)):
            return c
    return candidates[-1]


def BASALT_main_re_assembly(assembly_list, datasets, num_threads, lr_list, hifi_list,
                              hic_list, eb_list, ram, continue_mode, functional_module,
                              sensitive, refinement_paramter, max_ctn, min_cpn,
                              pwd, QC_software, output_folder):
    """
    Run the BASALT reassembly module. The QC backend (CheckM2 or

    CheckM) is selected at runtime via the ``QC_software`` parameter.
    All step modules are configured accordingly at function entry.
    """
    # Configure QC backend for all step modules in this phase.
    _set_s4_backend(QC_software)
    _set_s8_backend(QC_software)
    _set_s9_backend(QC_software)
    _set_s9p_backend(QC_software)
    _set_s10_backend(QC_software)

    ### Record the last accomplished step
    pwd=os.getcwd()

    #### Check existence of models
    user_dir = os.path.expanduser('~')
    # local_dir = f"{user_dir}/.cache/BASALT"
    BASALT_WEIGHT = os.environ.get("BASALT_WEIGHT")
    local_dir = BASALT_WEIGHT
    os.chdir(local_dir)
    model_list=glob(r'*_ensemble.csv')
    os.chdir(pwd)
    # print(model_list)
    if len(model_list) == 5:
        x=0
    else:
        print('BASALT models lacking. Start download the model')
        # Legacy standalone downloader replaced by basalt.ml.download.
        from basalt.ml.download import main as _download_models
        _download_models()

    #### Program start
    last_step=0
    if continue_mode == 'last':
        try:
            n=0
            for line in open('Basalt_checkpoint.txt', 'r'):
                n+=1

            n1=0
            for line in open('Basalt_checkpoint.txt', 'r'):
                n1+=1
                if n1 == n:
                    ls=str(line)[0]
                    try:
                        ls2=int(str(line)[1])
                        last_step=int(str(ls)+str(ls2))
                    except:
                        last_step=int(ls)
                    # last_step=int(str(line).replace('th','').replace('st','').replace('nd','').replace('rd','').split(' ')[0])
        except:
            f_cp_m=open('Basalt_checkpoint.txt', 'w')
            f_cp_m.close()
    else:
        # Reassembly consumes upstream artifacts produced by autobinning /
        # refinement. The global cleanup() would remove those, so a
        # downstream phase must NOT call it — only reset its own checkpoint.
        print('Resetting reassembly checkpoint (upstream outputs preserved)')
        f_cp_m=open('Basalt_checkpoint.txt', 'w')
        f_cp_m.close()

    print('BASALT started from step: '+str(last_step))
    try:
        fx=open('Basalt_log.txt','a')
        fx.write('BASALT started from step: '+str(last_step)+'\n')
        fx.close()
    except:
        fx=open('Basalt_log.txt','w')
        fx.write('BASALT started from step: '+str(last_step)+'\n')
        fx.close()
    fx.close()
    
    try:
        f=open('BASALT_command.txt','a')
        f.write('BASALT restarted from step: '+str(last_step)+'\n')
        f.close()
        fx=open('Basalt_log.txt','a')
        fx.write('BASALT restarted from step: '+str(last_step)+'\n')
        fx.close()
    except:
        f=open('BASALT_command.txt','w')
        f.write('BASALT started from 1st step'+'\n')
        f.write('Assemblies: '+str(str(assembly_list).replace('[','').replace(']','').replace(' ','').replace('\'',''))+'\n')
        f.write('Datasets: ')
        fx=open('Basalt_log.txt','a')
        fx.write('BASALT started from 1st step'+'\n')
        fx.write('Assemblies: '+str(str(assembly_list).replace('[','').replace(']','').replace(' ','').replace('\'',''))+'\n')
        fx.write('Datasets: ')
        for n_pair, ds_key in enumerate(datasets, start=1):
            pair = ','.join(datasets[ds_key])
            sep = '' if n_pair == 1 else '/'
            f.write(sep + pair)
            fx.write(sep + pair)
        f.write('\n')
        fx.write('\n')
        try:
            f.write('Long reads: '+str(str(lr_list).replace('[','').replace(']','').replace(' ','').replace('\'',''))+'\n')
            fx.write('Long reads: '+str(str(lr_list).replace('[','').replace(']','').replace(' ','').replace('\'',''))+'\n')
        except:
            f.write('Long reads: []'+'\n')
            fx.write('Long reads: []'+'\n')
        f.write('Num threads: '+str(num_threads)+'\n'+'Ram: '+str(ram)+'\n'+'Refinement mode: '+str(refinement_paramter)+'\n')
        f.write('Autobinning mode: '+str(sensitive)+'\n'+'Functional mode: '+str(functional_module)+'\n'+'Continue mode: '+str(continue_mode)+'\n')
        f.write('Max contamination: '+str(max_ctn)+'\n'+'Min completeness: '+str(min_cpn)+'\n')
        f.close()
        fx.write('Num threads: '+str(num_threads)+'\n'+'Ram: '+str(ram)+'\n'+'Refinement mode: '+str(refinement_paramter)+'\n')
        fx.write('Autobinning mode: '+str(sensitive)+'\n'+'Functional mode: '+str(functional_module)+'\n'+'Continue mode: '+str(continue_mode)+'\n')
        fx.write('Max contamination: '+str(max_ctn)+'\n'+'Min completeness: '+str(min_cpn)+'\n')
        fx.close()

    x = 0
    datasets2=copy.deepcopy(datasets)
    for item in datasets.keys():
        hz_list=datasets[item][0].split('.')
        if len(hz_list) >= 2:
            if hz_list[-1] == 'fq' or hz_list[-1] == 'fastq':
                x=1
            elif hz_list[-1] == 'zip':
                x=2
            elif hz_list[-1] == 'gz':
                if hz_list[-2] == 'tar':
                    x=3
                else:
                    x=4
        else:
            print('Input format error! Please check the input file.')
            print('BASALT supports the input (1) sequence files in  .gz, .zip, and .tar.gz; (2) and assemblies in .fa, .fna, .fasta.')
            fx=open('Basalt_log.txt','a')
            fx.write('Input format error! Please check the input file.'+'\n')
            fx.write('BASALT supports the input (1) sequence files in  .gz, .zip, and .tar.gz; (2) and assemblies in .fa, .fna, .fasta.'+'\n')
            fx.close()

    if x > 1:
        datasets={}
        for item in datasets2.keys():
            datasets[item]=[]
            if x == 2:
                f1_d=str(datasets2[item][0]).split('.zip')[0]
                f2_d=str(datasets2[item][1]).split('.zip')[0]
                if os.path.exists(pwd+'/PE_r1_'+str(f1_d)):
                    z=0
                else:
                    os.system('unzip '+str(datasets2[item][0]))

                if os.path.exists(pwd+'/PE_r2_'+str(f2_d)):
                    z=0
                else:
                    os.system('unzip '+str(datasets2[item][1]))

            elif x == 3:
                f1=str(datasets2[item][0])
                f2=str(datasets2[item][1])
                f1_d=str(f1).split('.tar.gz')[0]
                f2_d=str(f2).split('.tar.gz')[0]
                if os.path.exists(pwd+'/PE_r1_'+str(f1_d)):
                    z=0
                else:
                    os.system('tar -zxf '+str(f1))
                if os.path.exists(pwd+'/PE_r2_'+str(f2_d)):
                    z=0
                else:
                    os.system('tar -zxf '+str(f2))

            elif x == 4:
                f1=str(datasets2[item][0])
                f2=str(datasets2[item][1])
                f1_d=f1.split('.gz')[0]
                f2_d=f2.split('.gz')[0]
                if os.path.exists(pwd+'/PE_r1_'+str(f1_d)):
                    z=0
                else:
                    os.system('gunzip -c '+f1+' > '+f1_d)
                if os.path.exists(pwd+'/PE_r2_'+str(f2_d)):
                    z=0
                else:                    
                    os.system('gunzip -c '+f2+' > '+f2_d)

            if '.fq' not in f1_d and '.fastq' not in f1_d:
                f1_d=f1_d+'.fq'
            if '.fq' not in f2_d and '.fastq' not in f2_d:
                f2_d=f2_d+'.fq'
            datasets[item].append(f1_d)
            datasets[item].append(f2_d)

    x=0
    assembly_list2=copy.deepcopy(assembly_list)
    for item in assembly_list:
        hz_list=str(item).split('.')
        if len(hz_list) >= 2:
            if hz_list[-1] == 'fa' or hz_list[-1] == 'fasta' or hz_list[-1] == 'fna':
                x=1
            elif hz_list[-1] == 'zip':
                x=2
            elif hz_list[-1] == 'gz':
                if hz_list[-2] == 'tar':
                    x=3
                else:
                    x=4
        else:
            print('Input format error! Please check the input file.')
            print('BASALT supports the input (1) sequence files in  .gz, .zip, and .tar.gz; (2) and assemblies in .fa, .fna, .fasta.')
            fx=open('Basalt_log.txt','a')
            fx.write('Input format error! Please check the input file.'+'\n')
            fx.write('BASALT supports the input (1) sequence files in  .gz, .zip, and .tar.gz; (2) and assemblies in .fa, .fna, .fasta.'+'\n')
            fx.close()
    
    if x > 1:
        assembly_list=[]
        for item in assembly_list2:
            if x == 2:
                f_d=str(item).split('.zip')[0]
                if os.path.exists(pwd+'/'+str(f_d)):
                    z=0
                else:   
                    os.system('unzip '+str(item))
            
            elif x == 3:
                f_d=str(item).split('.tar.gz')[0]
                if os.path.exists(pwd+'/'+str(f_d)):
                    z=0
                else:
                    os.system('tar -zxf '+str(item))

            elif x == 4:
                f_d=str(item).split('.gz')[0]
                if os.path.exists(pwd+'/'+str(f_d)):
                    z=0
                else:
                    os.system('gunzip -c '+str(item)+' > '+str(f_d))

            if '.fa' not in f_d and '.fna' not in f_d and '.fasta' not in f_d:
                f_d=f_d+'.fa'
            assembly_list.append(f_d)

    ### Autobinner
    if functional_module == 'reassembly' or functional_module == 'all':
        if last_step < 8:
            if len(assembly_list) > 1:
                print('Starting OLC process')
                if len(lr_list) == 0 and len(hifi_list) == 0: ### only HTS datasets presented, perform OLC
                    # Resolve the binset fed into OLC, tolerating skipped upstream
                    # steps: on low-connectivity pure-SR input the contig-retrieve
                    # step (gated on PE/long-read connections in refinement) does
                    # not run, so *_filtrated_retrieved is never created. Pick the
                    # first existing non-empty binset in preference order and fall
                    # back to BestBinset_outlier_refined (always produced by S4),
                    # instead of handing OLC_main a missing folder (FileNotFoundError).
                    olc_candidates = [
                        'BestBinset_outlier_refined_filtrated_retrieved_retrieved',
                        'BestBinset_outlier_refined_filtrated_retrieved',
                        'BestBinset_outlier_refined',
                    ]
                    target_bin_folder = None
                    for _cand in olc_candidates:
                        if _folder_has_fasta(os.path.join(pwd, _cand)):
                            target_bin_folder = _cand
                            break
                    if target_bin_folder is None:
                        target_bin_folder = olc_candidates[-1]
                        print('No non-empty pre-OLC binset found; falling back to '+target_bin_folder)
                    else:
                        print('Selected pre-OLC binset: '+target_bin_folder)

                    print('Processing with bins in '+str(target_bin_folder))
                    bin_comparison_folder='BestBinset_comparison_files'
                    step='assemblies_OLC'
                    aligned_len_cutoff=500
                    similarity_cutoff=99
                    coverage_extension=95
                    if QC_software == 'checkm':
                        mod_bin_folder=''
                        OLC_main_checkm(target_bin_folder, step, bin_comparison_folder, aligned_len_cutoff, similarity_cutoff, coverage_extension, num_threads, ram, mod_bin_folder)
                    else:
                        mod_bin_folder=[]
                        OLC_main(target_bin_folder, step, bin_comparison_folder, aligned_len_cutoff, similarity_cutoff, coverage_extension, num_threads, mod_bin_folder)
                    f_cp_m=open('Basalt_checkpoint.txt', 'a')
                    f_cp_m.write('\n'+'8th contig OLC done!'+'\t'+target_bin_folder+'_OLC')
                    f_cp_m.close()
                else:
                    target_bin_folder=None
                    for line in open('Basalt_checkpoint.txt', 'r'):
                        if '7th' == str(line[0:3]):
                            target_bin_folder=str(line).strip().split('\t')[1].strip()
                    if target_bin_folder is None:
                        target_bin_folder=_resolve_pre_olc_bin_folder(pwd, lr_list)
                        print('No 7th-step entry in checkpoint; falling back to '+target_bin_folder)

                    print('Long-read presented. Skip OLC process.')
                    f_cp_m=open('Basalt_checkpoint.txt', 'a')
                    f_cp_m.write('\n'+'8th contig OLC did not perform!'+'\t'+target_bin_folder)
                    f_cp_m.close()
            else:
                # if len(lr_list) != 0: ### once lr_list presented, MAG polished present
                #     target_bin_folder='BestBinset_outlier_refined_filtrated_retrieved_MAGs_polished'
                # else: ### lr_list did not present. Negatelect polishing
                #     try:
                #         bin_n=0
                #         os.chdir(pwd+'/BestBinset_outlier_refined_filtrated_retrieved_retrieved')
                #         for root, dirs, files in os.walk(pwd+'/BestBinset_outlier_refined_filtrated_retrieved_retrieved'):
                #             for file in files:
                #                 hz=file.split('.')[-1]
                #                 if 'fa' in hz or 'fna' in hz:
                #                     bin_n+=1
                #         os.chdir(pwd)
                #         if bin_n != 0:
                #             target_bin_folder='BestBinset_outlier_refined_filtrated_retrieved_retrieved'
                #         else:
                #             target_bin_folder='BestBinset_outlier_refined_filtrated_retrieved'
                #     except:
                #         print('There is not bin in BestBinset_outlier_refined_filtrated_retrieved_retrieved. Try BestBinset_outlier_refined_filtrated_retrieved')
                #         target_bin_folder='BestBinset_outlier_refined_filtrated_retrieved'

                target_bin_folder=None
                for line in open('Basalt_checkpoint.txt', 'r'):
                    if '7th' == str(line[0:3]):
                        target_bin_folder=str(line).strip().split('\t')[1].strip()
                if target_bin_folder is None:
                    target_bin_folder=_resolve_pre_olc_bin_folder(pwd, lr_list)
                    print('No 7th-step entry in checkpoint; falling back to '+target_bin_folder)

                print('Long-read presented. Skip OLC process.')
                f_cp_m=open('Basalt_checkpoint.txt', 'a')
                f_cp_m.write('\n'+'8th contig OLC did not perform!'+'\t'+target_bin_folder)
                f_cp_m.close()

        if last_step < 9:
            if len(datasets) != 0: ### reassembly only when HTS presented
                print('Starting reassembly')
                for line in open('Basalt_checkpoint.txt', 'r'):
                    if '8th contig OLC' in line:
                        binset_folder=str(line).strip().split('\t')[1].strip()
                        print('Processing bins in '+str(binset_folder))

                datasets_list={}
                for ds in datasets.keys():
                    datasets_list[ds]=[]
                    datasets_list[ds].append('PE_r1_'+str(datasets[ds][0]))
                    datasets_list[ds].append('PE_r2_'+str(datasets[ds][1]))
                
                if len(lr_list) == 0:
                    hybri_reassembly='n'
                    re_assembly_main(binset_folder, datasets_list, lr_list, hybri_reassembly, ram, num_threads, sensitive=sensitive)
                    f_cp_m=open('Basalt_checkpoint.txt', 'a')
                    f_cp_m.write('\n'+'9th reassembly done.'+'\t'+str(binset_folder)+'_re-assembly')
                    f_cp_m.close()
                else:
                    polished_binset=binset_folder
                    sr_folder='BestBinset_sr_bins_seq'
                    lr_folder='BestBinset_long_read'
                    # polished_binset='BestBinset_outlier_refined_filtrated_retrieved_MAGs_polished'
                    # sr_folder='BestBinset_outlier_refined_filtrated_retrieved_polished_sr_bins_seq'
                    # lr_folder='BestBinset_outlier_refined_filtrated_retrieved_long_read'
                    hybrid_re_assembly_main(polished_binset, sr_folder, lr_folder, ram, num_threads, sensitive=sensitive)
                    os.system('rm -rf SPAdes_corrected_reads')
                    f_cp_m=open('Basalt_checkpoint.txt', 'a')
                    f_cp_m.write('\n'+'9th reassembly done.'+'\t'+str(polished_binset)+'_re-assembly')
                    f_cp_m.close()
            else:
                f_cp_m=open('Basalt_checkpoint.txt', 'a')
                f_cp_m.write('\n'+'9th Skipped the reassembly step.')
                f_cp_m.close()

        if last_step < 10:
            if len(datasets) != 0:
                print('Starting reassemblies OLC')
                for line in open('Basalt_checkpoint.txt', 'r'):
                    if '8th contig OLC' in line:
                        orig_binset=str(line).strip().split('\t')[1].strip()
                    elif '8th contig OLC did not perform' in line:
                        orig_binset=str(line).strip().split('\t')[1].strip()+'_mod'

                    if '9th reassembly done' in line:
                        target_bin_folder=str(line).strip().split('\t')[1].strip()
                        print('Processing bins in '+str(target_bin_folder))

                bin_comparison_or_reassembly_binset_folder=target_bin_folder+'_binset'
                # orig_binset=target_bin_folder.split('_re-assembly')[0]+'_mod'
                step='OLC_after_reassembly'
                aligned_len_cutoff=500
                similarity_cutoff=98
                coverage_extension=90

                reassembly_OLC_main(target_bin_folder, step, bin_comparison_or_reassembly_binset_folder, aligned_len_cutoff, similarity_cutoff, coverage_extension, num_threads, ram, orig_binset)
                f_cp_m=open('Basalt_checkpoint.txt', 'a')
                f_cp_m.write('\n'+'10th reassemblies OLC done.'+'\t'+target_bin_folder+'_OLC')
                f_cp_m.close()
            else:
                f_cp_m=open('Basalt_checkpoint.txt', 'a')
                f_cp_m.write('\n'+'10th Skipped reassemblies OLC.')
                f_cp_m.close()
        
        try:
            if last_step < 11:
                if len(lr_list) != 0 and len(datasets) != 0:
                    ### Long read present. Start polishing
                    print('Starting 2nd run polishing')
                    target_bin_folder='BestBinset_outlier_refined_filtrated_retrieved_MAGs_polished_re-assembly_OLC'
                    remained_seq_list=glob(r'Remained_seq1.fq_*')
                    assembly_mo_list = []
                    for line in open('Assembly_mo_list.txt','r'):
                        assembly_mo_list.append(str(line).strip())

                    i=0
                    for item in remained_seq_list:
                        num=int(item.split('.fq_')[1])
                        if num > i:
                            i = num
                    os.system('mv Remained_seq1.fq_'+str(num)+' 2nd_Remained_seq1.fq')
                    os.system('mv Remained_seq2.fq_'+str(num)+' 2nd_Remained_seq2.fq')
                    unmapped_datasets={'1':['2nd_Remained_seq1.fq','2nd_Remained_seq2.fq']}
                    # unmapped_datasets={'1':['Remained_seq1.fq','Remained_seq2.fq']}


                    polishing_main(target_bin_folder, unmapped_datasets, assembly_mo_list, lr_list, 2, 'bwa', 3, ram, num_threads)
                    # polishing_main(target_bin_folder, unmapped_datasets, assembly_list, lr_list, 2, 'bwa', 3, ram, num_threads)
                    remained_seq_list=glob(r'Remained_seq1.fq_*')
                    i=0
                    for item in remained_seq_list:
                        num=int(item.split('.fq_')[1])
                        if num > i:
                            i = num
                    
                    for i in range(1,num+1):
                        os.system('rm Remained_seq1.fq_'+str(i)+' '+'Remained_seq2.fq_'+str(i))
                    os.system('mkdir Remained_seq')
                    os.system('mv 2nd_Remained_seq1.fq '+pwd+'/Remained_seq/Remained_seq1.fq')
                    os.system('mv 2nd_Remained_seq2.fq '+pwd+'/Remained_seq/Remained_seq2.fq')
                    os.system('tar -zcvf Remained_seq.tar.gz Remained_seq')
                    os.system('rm -rf Remained_seq')
                    f_cp_m=open('Basalt_checkpoint.txt', 'a')
                    f_cp_m.write('\n'+'11th 2nd Polishing done!')
                    f_cp_m.close()
                else:
                    f_cp_m=open('Basalt_checkpoint.txt', 'a')
                    f_cp_m.write('\n'+'11th did not perform polish')
                    f_cp_m.close()
        except:
            print('Did not perform the 2nd polish. Skip the step.')
            f_cp_m=open('Basalt_checkpoint.txt', 'a')
            f_cp_m.write('\n'+'11th 2nd polishing did not perform')
            f_cp_m.close()

        if last_step < 12:
            ### Final de-replication: S4 function
            for line in open('Basalt_checkpoint.txt', 'r'):
                try:
                    if '10th reassemblies OLC done' in line:
                        final_folder=str(line).strip().split('\t')[1].strip()
                except:
                    xyz=0

                try:
                    if '11th 2nd Polishing done' in line:
                        final_folder='BestBinset_outlier_refined_filtrated_retrieved_MAGs_polished_re-assembly_OLC_MAGs_polished'
                except:
                    xyz=0
                
                try:
                    if '11th 2nd did not perform' in line:
                        final_folder='BestBinset_outlier_refined_filtrated_retrieved_MAGs_polished_re-assembly_OLC'
                except:
                    xyz=0
                
                try:
                    if len(lr_list) != 0 and len(datasets_list) == 0: ### only long reads, excluding hifi
                        final_folder='BestBinset_outlier_refined_MAGs_polished'
                except:
                    xyz=0
                
                try:
                    if len(hifi_list) != 0 and len(datasets_list) == 0 and len(lr_list) == 0: ### only hifi
                        final_folder='BestBinset_outlier_refined_filtrated_retrieved'
                except:
                    xyz=0

            if len(assembly_list) != 1:                
                coverage_matrix_list, bestbinset_list = [], []
                # for line in open('Bestbinset_list.txt','r'):
                #     bestbinset_list.append(str(line).strip())
            
                # for line in open('Coverage_matrix_list.txt','r'):
                #     coverage_matrix_list.append(str(line).strip())

                # datasets_fq={}
                # for item in datasets.keys():
                #     datasets_fq[item]=[]
                #     datasets_fq[item].append('PE_r1_'+str(datasets[item][0]))
                #     datasets_fq[item].append('PE_r2_'+str(datasets[item][1]))

                # multiple_assembly_comparator_main(drep_list, bestbinset_list, coverage_matrix_list, datasets_fq, 'final_drep', num_threads)
                try:
                    final_binset_comparator(final_folder, coverage_matrix_list, datasets, num_threads, pwd, 'final_drep')
                except:
                    print('Final drep did not perform')
            #     # final_binset_comparator(final_folder, coverage_matrix_list, datasets_fq, num_threads, pwd, 'final_drep')
            f_cp_m=open('Basalt_checkpoint.txt', 'a')
            f_cp_m.write('\n'+'12th final de-replication done!')
            f_cp_m.close()

        ### Final binset
        os.system('mv '+str(final_folder)+' '+str(output_folder))

        rename={}
        os.chdir(pwd+'/'+str(output_folder))
        rename_file=open('Rename.txt', 'w')
        for root, dirs, files in os.walk(pwd+'/'+str(output_folder)):
            for file in files:
                hz=file.split('.')[-1]
                if 'fa' in hz or 'fna' in hz or 'fasta' in hz:
                    if '_' in file:
                        nl=file.split('.')
                        nl.remove(nl[-1])
                        full_name='.'.join(nl)
                        bin_name=file.split('_')[0].split('.')[0]
                        os.system('mv '+file+' '+bin_name+'.fa')
                        rename[full_name]=bin_name
                        rename_file.write(file+' '+bin_name+'.fa'+'\n')
        rename_file.close()

        x=0
        for root, dirs, files in os.walk(pwd+'/'+str(output_folder)):
            for file in files:
                if 'quality_report.tsv' in file:
                    x=1
                    f=open(file+'2','w')
                    xyz=0
                    for line in open(file):
                        xyz+=1
                        if xyz >= 2:
                            n=str(line).strip().split('\t')[0]
                            # c=str(line).strip().split('\t')[1]
                            if n in rename.keys():
                                f.write(line.replace(n, str(rename[n])))
                            else:
                                f.write(line)
                        else:
                            f.write(line)
                    f.close()
                    
                    os.system('mv '+file+'2 '+file)
        
        os.system('rm *quality_report_o.tsv')
        os.chdir(pwd)

        if x == 0:
            os.system('checkm2 predict -t '+str(num_threads)+' -i '+str(output_folder)+' -x fa -o Final_bestbinset_checkm')
            os.system('mv '+pwd+'/Final_bestbinset_checkm/quality_report.tsv '+pwd+'/'+str(output_folder)+'/Final_bestbinset_quality_report.tsv')
        else:
            print('Final quality report existed')

        f_cp_m=open('Basalt_checkpoint.txt', 'a')
        f_cp_m.write('\n'+'BASALT done!')
        f_cp_m.close()
    print('BASALT main program accomplished!')
    print('BASALT will continue to cleanup or compress all the temp files. The results could be found in folder '+repr(str(output_folder))+'. Please wait for a little bit longer')

    ### Cleanup

if __name__ == '__main__':
    assembly_list=['8_medium_S001_SPAdes_scaffolds.fasta','10_medium_cat_SPAdes_scaffolds.fasta']
    datasets={'1':['RM2_S001_insert_270_mate1.fq','RM2_S001_insert_270_mate2.fq']}
    lr_list=[]
    hic_list=[]
    hifi_list=[]
    eb_list=[]
    num_threads=30
    ram=120
    pwd=os.getcwd()
    refinement_paramter='deep'
    sensitive='sensitive'
    functional_module='reassembly'
    continue_mode='last'
    max_ctn, min_cpn=20, 35
    QC_software='checkm2'
    output_folder='Final_binset'
    BASALT_main_re_assembly(assembly_list, datasets, num_threads, lr_list, hifi_list, hic_list, eb_list, ram, continue_mode, functional_module, sensitive, refinement_paramter, max_ctn, min_cpn, pwd, QC_software, output_folder)
