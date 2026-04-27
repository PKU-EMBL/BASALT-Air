#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Refinement entry point for the BASALT pipeline.

This module runs DL-based outlier removal (S5) and within-group contig
retrieval (S6 / S7 / S7lr), records progress in ``Basalt_checkpoint.txt``
to support resumable execution.

Supports both CheckM2 and CheckM via the ``QC_software`` parameter.
"""

import time
import sys
import os
from Bio import SeqIO
from basalt.steps.s4_multiple_assembly_comparator import *
from basalt.steps.s4_multiple_assembly_comparator import set_qc_backend as _set_s4_backend
from basalt.steps.s5_outlier_remover_dl import *
from basalt.steps.s5_outlier_remover_dl import set_qc_backend as _set_s5_backend
from basalt.steps.s6_retrieve_contigs_from_pe_contigs import *
from basalt.steps.s6_retrieve_contigs_from_pe_contigs import set_qc_backend as _set_s6_backend
from basalt.steps.s7_contigs_retrieve_within_group import *
from basalt.steps.s7_contigs_retrieve_within_group import set_qc_backend as _set_s7_backend
from basalt.steps.s7lr_finding_sr_contigs import *
from basalt.steps.s7lr_finding_sr_contigs import set_qc_backend as _set_s7lr_backend
from glob import glob
from basalt.core.cleanup import *


def BASALT_main_refinement(assembly_list, datasets, num_threads, lr_list, hifi_list,
                             hic_list, eb_list, ram, continue_mode, functional_module,
                             sensitive, refinement_paramter, max_ctn, min_cpn,
                             pwd, QC_software, output_folder):
    """
    Run the BASALT refinement module. The QC backend (CheckM2 or

    CheckM) is selected at runtime via the ``QC_software`` parameter.
    All step modules are configured accordingly at function entry.
    """
    # Configure QC backend for all step modules in this phase.
    _set_s4_backend(QC_software)
    _set_s5_backend(QC_software)
    _set_s6_backend(QC_software)
    _set_s7_backend(QC_software)
    _set_s7lr_backend(QC_software)

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
        # os.system('python BASALT_models_download.py')
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
        print('Start a new project')
        cleanup(assembly_list)
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
        datasets_list=datasets.split(',')
        x=0
        for ds in datasets_list:
            x+=1
            ds1=str(ds).split('[')[1].replace('\'','').replace(']','').replace(' ','')
            if x == 1:
                f.write(str(ds1))
                fx.write(str(ds1))
            else:
                f.write('/'+str(ds1))
                fx.write('/'+str(ds1))
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
            print('BASALT supports the input (1) sequence files in  .gz, .zip, and .tar.gz; (2) and assemlies in .fa, .fna, .fasta.')
            fx=open('Basalt_log.txt','a')
            fx.write('Input format error! Please check the input file.'+'\n')
            fx.write('BASALT supports the input (1) sequence files in  .gz, .zip, and .tar.gz; (2) and assemlies in .fa, .fna, .fasta.'+'\n')
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
        hz_list=assembly_list[0].split('.')
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
            print('BASALT supports the input (1) sequence files in  .gz, .zip, and .tar.gz; (2) and assemlies in .fa, .fna, .fasta.')
            fx=open('Basalt_log.txt','a')
            fx.write('Input format error! Please check the input file.'+'\n')
            fx.write('BASALT supports the input (1) sequence files in  .gz, .zip, and .tar.gz; (2) and assemlies in .fa, .fna, .fasta.'+'\n')
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
    if functional_module == 'refinement' or functional_module == 'all':
        if last_step < 4:
            print('Starting outlier removal process')
            coverage_matrix_list, connections_list, assembly_mo_list, bestbinset_list = [], [], [], []
            for line in open('Coverage_matrix_list.txt','r'):
                coverage_matrix_list.append(str(line).strip())
            for line in open('Assembly_mo_list.txt','r'):
                assembly_mo_list.append(str(line).strip())
            for line in open('Bestbinset_list.txt','r'):
                bestbinset_list.append(str(line).strip())
            
            for item in assembly_list:
                connections_list.append('condense_connections_'+item+'.txt')
            if len(bestbinset_list) == 1:
                print('Copying '+str(bestbinset_list[0])+' to BestBinset')
                os.system('cp -r '+str(bestbinset_list[0])+' BestBinset')
            # contig_outlier_remover_main('BestBinset', coverage_matrix_list, connections_list, num_threads, ram)
            # outlier_remover_main('BestBinset', coverage_matrix_list, datasets, assembly_mo_list, pwd, num_threads)
            # lr_list2=copy.deepcopy(lr_list)
            # for item2 in hifi_list:
            #     lr_list2.append(item2)            

            outlier_remover_main('BestBinset', coverage_matrix_list, datasets, assembly_mo_list, pwd, num_threads, lr=lr_list, hifi_list=hifi_list)
            f_cp_m=open('Basalt_checkpoint.txt', 'a')
            f_cp_m.write('\n'+'4th outlier removal done!')
            f_cp_m.close()

        if last_step < 5:
            print('Starting contig retrival process')
            best_binset_from_multi_assemblies='BestBinset_outlier_refined'
            outlier_remover_folder='BestBinset_outlier_refined'
            
            coverage_matrix_list, connections_list, assembly_mo_list = [], [], []
            for line in open('Coverage_matrix_list.txt','r'):
                coverage_matrix_list.append(str(line).strip())
            for line in open('Assembly_mo_list.txt','r'):
                assembly_mo_list.append(str(line).strip())
            for line in open('Assembly_MoDict.txt','r'):
                item=str(line).strip().split('\t')[0].strip()
                connections_list.append('condense_connections_'+item+'.txt')

            xxx=0
            try:
                
                for line in open('Long_reads_connecting_contigs.txt','r'):
                    xxx+=1
            except:
                xxx=0
            
            yyy=0
            for item in connections_list:
                xyz=0
                for line in open(item,'r'):
                    xyz+=1
                if xyz > 1:
                    yyy+=1

            if xxx >= 1:
                lr_connection_list='Long_reads_connecting_contigs.txt'
            else:
                lr_connection_list=''
            
            if xxx != 0 or yyy != 0:
                if sensitive == 'more-sensitive':
                    refinement_paramter = 'deep'
                # else:
                #     refinement_paramter = 'quick'
                Contig_recruiter_main(best_binset_from_multi_assemblies, outlier_remover_folder, num_threads, continue_mode, min_cpn, max_ctn, assembly_mo_list, connections_list, coverage_matrix_list, refinement_paramter, pwd, lr_connection_list=lr_connection_list)
                # Contig_recruiter_main(best_binset_from_multi_assemblies, outlier_remover_folder, num_threads, parameter, cpn_cutoff, ctn_cutoff, assemblies_list, PE_connections_list, lr_connection_list, coverage_matrix_list, refinement_mode, pwd)
                try:
                    os.chdir(str(best_binset_from_multi_assemblies)+'_retrieved_checkm')
                    n1=0
                    for line in open('quality_report.tsv','r'):
                        n1+=1
                        if n1 == 2:
                            break
                    os.chdir(str(pwd))
                
                    f_cp_m=open('Basalt_checkpoint.txt', 'a')
                    f_cp_m.write('\n'+'5th contig retrieve done!')
                    f_cp_m.close()
                except:
                    f_cp_m=open('Basalt_checkpoint.txt', 'a')
                    f_cp_m.write('\n'+'5th contig retrieve did not perform!')
                    f_cp_m.close()
                
            elif xxx == 0 or yyy == 0:
                f_cp_m=open('Basalt_checkpoint.txt', 'a')
                f_cp_m.write('\n'+'5th contig retrieve did not perform!')
                f_cp_m.close()

            os.system('cp -r BestBinset_outlier_refined_filtrated_retrieved BestBinset_outlier_refined_filtrated_retrieved_backup')

        if last_step < 6:
            ### Second de-replication: S4 function
            if len(assembly_list) != 1:
                print('Secondary de-repplication start')
                #drep_list='BestBinset_outlier_refined_filtrated_retrieved'

                coverage_matrix_list, bestbinset_list = [], []
                for line in open('Bestbinset_list.txt','r'):
                    bestbinset_list.append(str(line).strip())
            
                for line in open('Coverage_matrix_list.txt','r'):
                    coverage_matrix_list.append(str(line).strip())

                for line in open('Basalt_checkpoint.txt','r'):
                    if line[0] == '5':
                        if '5th contig retrieve done!' in line:
                            d_folder='BestBinset_outlier_refined_filtrated_retrieved'
                        elif '5th contig retrieve did not perform!' in line:
                            d_folder='BestBinset_outlier_refined'

                datasets_fq={}
                for item in datasets.keys():
                    datasets_fq[item]=[]
                    datasets_fq[item].append('PE_r1_'+str(datasets[item][0]))
                    datasets_fq[item].append('PE_r2_'+str(datasets[item][1]))

                # multiple_assembly_comparator_main(drep_list, bestbinset_list, coverage_matrix_list, datasets_fq, 'second_drep', num_threads)
                # final_binset_comparator('BestBinset_outlier_refined_filtrated_retrieved', coverage_matrix_list, datasets_fq, num_threads, pwd, 'second_drep')
                final_binset_comparator(d_folder, coverage_matrix_list, datasets_fq, num_threads, pwd, 'second_drep')
            f_cp_m=open('Basalt_checkpoint.txt', 'a')
            f_cp_m.write('\n'+'6th secondary de-repplication done!')
            f_cp_m.close()

        if last_step < 7:
            if len(lr_list) == 0 and len(hifi_list) == 0: ### only HTS datasets presented
                # print('Starting contig retrival within group')
                # Contig_recruiter_main('BestBinset_outlier_refined', num_threads, parameter, 35, 20)
                best_binset_after_contig_retrieve='BestBinset_outlier_refined_filtrated_retrieved'
                outlier_remover_folder='BestBinset_outlier_refined'
                cpn_cutoff=35 ### The minimun completeness to keep
                ctn_cutoff=5 ### The maximal contaminaition to keep

                coverage_matrix_list, connections_list, assembly_mo_list = [], [], []
                for line in open('Coverage_matrix_list.txt','r'):
                    coverage_matrix_list.append(str(line).strip())
                for line in open('Assembly_mo_list.txt','r'):
                    assembly_mo_list.append(str(line).strip())
                for line in open('Assembly_MoDict.txt','r'):
                    item=str(line).strip().split('\t')[0].strip()
                    connections_list.append('condense_connections_'+item+'.txt')

                # Contig_retrieve_within_group_main(best_binset_after_contig_retrieve, outlier_remover_folder, num_threads, continue_mode, cpn_cutoff, ctn_cutoff, assembly_mo_list, connections_list, coverage_matrix_list)

                f_cp_m=open('Basalt_checkpoint.txt', 'a')
                f_cp_m.write('\n'+'7th Skip contig retrieve within group!'+'\t'+best_binset_after_contig_retrieve)
                # f_cp_m.write('\n'+'7th contig retrieve within group done!'+'\t'+best_binset_after_contig_retrieve+'_retrieved')
                f_cp_m.close()
                 
            if len(hifi_list) != 0 and len(lr_list) == 0 and len(datasets) == 0: ### only hifi data presented, skip polishing
                best_binset_after_contig_retrieve='BestBinset_outlier_refined_filtrated_retrieved'
                outlier_remover_folder='BestBinset_outlier_refined'
                cpn_cutoff=35 ### The minimun completeness to keep
                ctn_cutoff=5 ### The maximal contaminaition to keep

                coverage_matrix_list, connections_list, assembly_mo_list = [], [], []
                for line in open('Coverage_matrix_list.txt','r'):
                    coverage_matrix_list.append(str(line).strip())
                for line in open('Assembly_mo_list.txt','r'):
                    assembly_mo_list.append(str(line).strip())
                for line in open('Assembly_MoDict.txt','r'):
                    item=str(line).strip().split('\t')[0].strip()
                    connections_list.append('condense_connections_'+item+'.txt')

                # Contig_retrieve_within_group_main(best_binset_after_contig_retrieve, outlier_remover_folder, num_threads, continue_mode, cpn_cutoff, ctn_cutoff, assembly_mo_list, connections_list, coverage_matrix_list)

                f_cp_m=open('Basalt_checkpoint.txt', 'a')
                f_cp_m.write('\n'+'7th Skip contig retrieve within group!'+'\t'+best_binset_after_contig_retrieve)
                # f_cp_m.write('\n'+'7th contig retrieve within group done!'+'\t'+best_binset_after_contig_retrieve+'_retrieve')
                f_cp_m.close()

            if len(lr_list) != 0:
                ### Long read present. Start polishing
                print('Starting 1st run polishing')
                # if len(datasets) != 0: 

                for line in open('Basalt_checkpoint.txt','r'):
                    if line[0] == '5':
                        if '5th contig retrieve done!' in line:
                            best_binset_after_contig_retrieve='BestBinset_outlier_refined_filtrated_retrieved'
                        elif '5th contig retrieve did not perform!' in line:
                            best_binset_after_contig_retrieve='BestBinset_outlier_refined'
                    
                assembly_mo_list = []
                for line in open('Assembly_mo_list.txt','r'):
                    assembly_mo_list.append(str(line).strip())

                if functional_module == 'all':
                    # best_binset_after_contig_retrieve='BestBinset_outlier_refined_filtrated_retrieved'
                    polishing_main(best_binset_after_contig_retrieve, datasets, assembly_mo_list, lr_list, 1, 'bw2', 3, ram, num_threads)
                    # polishing_main(best_binset_after_contig_retrieve, datasets, assembly_list, lr_list, 1, 'bw2', 3, ram, num_threads)
                    f_cp_m=open('Basalt_checkpoint.txt', 'a')
                    f_cp_m.write('\n'+'7th Polishing done!'+'\t'+best_binset_after_contig_retrieve+'_MAGs_polished')
                    f_cp_m.close()
                elif functional_module == 'refinement':
                    # best_binset_after_contig_retrieve='BestBinset_outlier_refined_filtrated_retrieved'
                    polishing_main(best_binset_after_contig_retrieve, datasets, assembly_mo_list, lr_list, 1, 'bwa', 3, ram, num_threads)
                    # polishing_main(best_binset_after_contig_retrieve, datasets, assembly_list, lr_list, 1, 'bwa', 3, ram, num_threads)
                    f_cp_m=open('Basalt_checkpoint.txt', 'a')
                    f_cp_m.write('\n'+'7th Polishing done!'+'\t'+best_binset_after_contig_retrieve+'_MAGs_polished')
                    f_cp_m.close()


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
    functional_module='refinement'
    continue_mode='last'
    max_ctn, min_cpn=20, 35
    QC_software='checkm2'
    output_folder='Final_binset'
    BASALT_main_refinement(assembly_list, datasets, num_threads, lr_list, hifi_list, hic_list, eb_list, ram, continue_mode, functional_module, sensitive, refinement_paramter, max_ctn, min_cpn, pwd, QC_software, output_folder)
