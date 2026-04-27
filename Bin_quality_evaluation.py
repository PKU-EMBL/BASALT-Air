#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#coding=utf-8

from Bio import SeqIO
import os, copy, math
import os, threading
from multiprocessing import Pool

def cat_contigs(bin_folder, pwd):
    bin_length, bin_contigs_n, bin_contigs_5k, bin_avg_length, bin_n50_length, bin_contigs_position_list={}, {}, {}, {}, {}, {}
    try:
        bin_folder_name=bin_folder.replace(' ','_')
    except:
        bin_folder_name=bin_folder

    contig_file=open('Total_contigs_'+bin_folder_name+'.fa','w')
    bin_summary_file=open(bin_folder_name+'_basic_information_of_bins.txt','w')
    bin_summary_file.write('Bin'+'\t'+'Bin total length (bp)'+'\t'+'NO. contigs'+'\t'+'NO. contigs with than 5kbp'+'\t'+'Bin avg. length (bp)'+'\t'+'Bin N50 contig\'s length (bp)'+'\n')
    os.chdir(pwd+'/'+bin_folder)
    for root, dirs, files in os.walk(pwd+'/'+bin_folder):
        for file in files:
            hz=str(file).split('.')[-1]
            if 'fa' in hz or 'fna' in hz:
                bin_contig_rank=[]
                bin_length[file]=0
                bin_contigs_n[file]=0
                bin_contigs_position_list[file], n, n1 = {}, 0, 0
                for record in SeqIO.parse(file, 'fasta'):
                    n+=1
                    bin_length[file]+=len(record.seq)
                    if len(record.seq) >= 5000:
                        n1+=1
                    bin_contig_rank.append(len(record.seq))
                    bin_contigs_position_list[file][record.id]={}
                    contig_file.write('>'+str(file)+'||'+str(record.id)+'\n'+str(record.seq)+'\n')
                bin_contig_rank.sort()
                bin_avg_length[file]=int(bin_length[file])/n
                bin_contigs_n[file], bin_contigs_5k[file]=n, n1

                n50_contig_index=math.ceil(float(n)/2)

                x=0
                for item in bin_contig_rank:
                    x+=1
                    if x == n50_contig_index:
                        bin_n50_length[file]=str(n50_contig_index)+':'+str(item)
    os.chdir(pwd)
    
    for bins in bin_length.keys():
        bin_summary_file.write(str(bins)+'\t'+str(bin_length[bins])+'\t'+str(bin_contigs_n[bins])+'\t'+str(bin_contigs_5k[bins])+'\t'+str(bin_avg_length[bins])+'\t'+str(bin_n50_length[bins])+'\n')
    contig_file.close()
    bin_summary_file.close()
    return 'Total_contigs_'+bin_folder_name+'.fa', bin_length

def alignment_len_calc(blast_output, bin_length):
    vs_bins_contigs, aligned_length={}, {}
    standard_bin=str(blast_output).split('_vs_')[0].split('Filtrated_')[1]
    test_bin=str(blast_output).split('_vs_')[1].split('.txt')[0]
    vs_bins_contigs[standard_bin]={}
    vs_bins_contigs[test_bin]={}

    for line in open(blast_output,'r'):
        simi=float(str(line).strip().split('\t')[2])
        length=float(str(line).strip().split('\t')[3])
        if float(simi) >= 99 and int(length) >= 100:
            contig1=str(line).strip().split('\t')[0].split('||')[1]  
            query_start=int(str(line).strip().split('\t')[6])
            query_end=int(str(line).strip().split('\t')[7])
            # standard_bin=str(line).strip().split('\t')[0].split('||')[0]       
            contig2=str(line).strip().split('\t')[1].split('||')[1] 
            subject_start=int(str(line).strip().split('\t')[8])
            subject_end=int(str(line).strip().split('\t')[9])
            # test_bin=str(line).strip().split('\t')[1].split('||')[0] 
            # vs_bins2=standard_bin+'\t'+test_bin
            # if vs_bins == vs_bins2:
            if contig1 not in vs_bins_contigs[standard_bin].keys():
                vs_bins_contigs[standard_bin][contig1]={}
    
            if contig2 not in vs_bins_contigs[test_bin].keys():
                vs_bins_contigs[test_bin][contig2]={}

            for i in range(query_start, query_end+1):
                vs_bins_contigs[standard_bin][contig1][i]=1
        
            if subject_end > subject_start:
                for i in range(subject_start, subject_end+1):
                    vs_bins_contigs[test_bin][contig2][i]=1
            else:
                for i in range(subject_end, subject_start+1):
                    vs_bins_contigs[test_bin][contig2][i]=1
    
    for bins in vs_bins_contigs.keys():
        aligned_length[bins]=0
        for contig in vs_bins_contigs[bins].keys():
            for position in vs_bins_contigs[bins][contig].keys():
                aligned_length[bins]+=vs_bins_contigs[bins][contig][position]
    # print(aligned_length)
    # print(str(bin_length[standard_bin]))
    output='Summary_'+str(standard_bin)+'_vs_'+str(test_bin)+'.txt'
    f=open(output,'w')
    aligned_ratio_standard_bin=100*float(aligned_length[standard_bin])/float(bin_length[standard_bin])
    aligned_ratio_test_bin=100*float(aligned_length[test_bin])/float(bin_length[test_bin])
    contamination=100*(float(bin_length[test_bin])-float(aligned_length[test_bin]))/float(bin_length[standard_bin])
    CPN_CTN=float(aligned_ratio_standard_bin)-float(contamination)
    CPN_CTN3=float(aligned_ratio_standard_bin)-3*float(contamination)
    CPN_CTN5=float(aligned_ratio_standard_bin)-5*float(contamination)
    # print(str(str(vs_bin)+'\t'+str(bin_length[standard_bin])+'\t'+str(bin_length[test_bin])))
    f.write(str(standard_bin)+'\t'+str(test_bin)+'\t'+str(bin_length[standard_bin])+'\t'+str(bin_length[test_bin]))
    f.write('\t'+str(aligned_length[standard_bin])+'\t'+str(aligned_length[test_bin])+'\t'+str(aligned_ratio_standard_bin)+'\t'+str(aligned_ratio_test_bin))
    f.write('\t'+str(contamination)+'\t'+str(CPN_CTN)+'\t'+str(CPN_CTN3)+'\t'+str(CPN_CTN5)+'\n')
    f.close()
    os.system('rm '+blast_output)

def ORF_predict(contig_file, pwd, standard_bin_orfs_folder):
    try:
        print('Predicting ORFs of '+str(contig_file))
        os.system('prodigal -i '+str(contig_file)+' -d '+str(contig_file)+'.orfs.fna -m -p meta -q')
        os.system('mv '+str(contig_file)+'.orfs.fna '+pwd+'/'+standard_bin_orfs_folder)
    except:
        print('ORFs prediction error! Please make sure prodigal is installed in your system')

def batch_makeblastdb(gs_bin):
    os.system('makeblastdb -in '+str(gs_bin)+'.orfs.fna -dbtype nucl')

def batch_blast(gs_bin, test_bin):
    os.system('blastn -db '+str(gs_bin)+'.orfs.fna -query '+str(test_bin)+'.orfs.fna -evalue 1e-20 -outfmt 6 -num_threads 1 -out '+str(gs_bin)+'_vs_'+(test_bin)+'.txt')

def Contigs_aligner(standard_bin_contigs, test_bin_contigs, bin_length, standard_bin_folder, test_bin_folder, num_threads, pwd):
    print('Using BLAST to align', test_bin_contigs, 'to', standard_bin_contigs)
    blast_output=str(test_bin_contigs)+'_vs_'+str(standard_bin_contigs)+'.txt'
    # os.system('makeblastdb -in '+str(test_bin_contigs)+' -dbtype nucl -hash_index -parse_seqids')
    os.system('makeblastdb -in '+str(test_bin_contigs)+' -dbtype nucl')
    os.system('blastn -query '+str(standard_bin_contigs)+' -db '+str(test_bin_contigs)+' -evalue 1e-20 -outfmt 6 -num_threads '+str(num_threads)+' -out '+str(blast_output))

    # blast_output2=open('Filtrated_'+str(test_bin_contigs)+'_vs_'+str(standard_bin_contigs)+'.txt','w')
    bin_paired_aligned_length, bin_paired_aligned_similarity, bin_paired_aligned_contamination_ratio={}, {}, {}
    vs_bins_contigs, output_files, vs_bin_job={}, {}, {}
    for line in open(blast_output,'r'):
        simi=float(str(line).strip().split('\t')[2])
        length=float(str(line).strip().split('\t')[3])
        if float(simi) >= 99 and int(length) >= 100:
            standard_bin=str(line).strip().split('\t')[0].split('||')[0]       
            test_bin=str(line).strip().split('\t')[1].split('||')[0] 
            vs_bins=standard_bin+'\t'+test_bin
            if vs_bins not in vs_bins_contigs.keys():
                f=open('Filtrated_'+standard_bin+'_vs_'+test_bin+'.txt','w')
                f.write(line)
                f.close()
                vs_bins_contigs[vs_bins]={}
                vs_bin_job[vs_bins]='Filtrated_'+standard_bin+'_vs_'+test_bin+'.txt'
                output_files['Summary_'+standard_bin+'_vs_'+test_bin+'.txt']=0
                vs_bins_contigs[vs_bins][standard_bin]={}
                vs_bins_contigs[vs_bins][test_bin]={}
                # aligned_length[vs_bins]={}
            else:
                f=open('Filtrated_'+standard_bin+'_vs_'+test_bin+'.txt','a')
                f.write(line)
                f.close()

    pool=Pool(processes=num_threads)
    x = 0
    for vs_bins in vs_bin_job.keys():
        x+=1
        blast_output=str(vs_bin_job[vs_bins])
        print('Processing '+str(blast_output))
        pool.apply_async(alignment_len_calc, args=(blast_output, bin_length))
    pool.close()
    pool.join()

    test_bin_non_redudant, test_bin_non_redudant2={}, {}
    f=open(str(test_bin_contigs)+'_bin_similarity_contamination_with_standard.txt','w')
    f.write('Standard bin'+'\t'+'Test bin'+'\t'+'Standard bin length(bp)'+'\t'+'Test bin aligned length(bp)'+'\t'+'Standard bin aligned length(bp)'+'\t'+'Test bin aligned length(bp)'+'\t'+'Standard bin aligned ratio(%)'+'\t'+'Test bin aligned ratio(%)'+'\t'+'Contaminated ratio(%)'+'\t'+'CPN-CTN(%)'+'\t'+'CPN-3*CTN(%)'+'\t'+'CPN-5*CTN(%)'+'\n')
    for output in output_files.keys():
        for line in open(output,'r'):
            f.write(str(line))
            gs_bin=str(line).strip().split('\t')[0]
            test_bin=str(line).strip().split('\t')[1]
            q_value=float(str(line).strip().split('\t')[-1])
            if test_bin in test_bin_non_redudant.keys():
                if q_value > test_bin_non_redudant[test_bin]:
                    test_bin_non_redudant[test_bin]=q_value
                    test_bin_non_redudant2[test_bin]=str(line)
            else:
                test_bin_non_redudant[test_bin]=q_value
                test_bin_non_redudant2[test_bin]=str(line)
        os.system('rm '+output)
    f.close()

    gs_redundant_bin={}
    f=open(test_bin_contigs+'_test_non-redundant_bin_similarity_contamination_with_standard.txt','w')
    f.write('Standard bin'+'\t'+'Test bin'+'\t'+'Standard bin length(bp)'+'\t'+'Test bin aligned length(bp)'+'\t'+'Standard bin aligned length(bp)'+'\t'+'Test bin aligned length(bp)'+'\t'+'Standard bin aligned ratio(%)'+'\t'+'Test bin aligned ratio(%)'+'\t'+'Contaminated ratio(%)'+'\t'+'CPN-CTN(%)'+'\t'+'CPN-3*CTN(%)'+'\t'+'CPN-5*CTN(%)'+'\n')
    for test_bin in test_bin_non_redudant2.keys():
        line=test_bin_non_redudant2[test_bin]
        gs_bin=line.strip().split('\t')[0]
        test_bin=line.strip().split('\t')[1]
        gs_redundant_bin[gs_bin]=0
        f.write(line)
    f.close()

    gs_redundant_bin2, bin_with_microdiversity, bin_with_microdiversity2={}, {},{}
    for test_bin in test_bin_non_redudant2.keys():
        line=test_bin_non_redudant2[test_bin]
        gs_bin=line.strip().split('\t')[0]
        q_value=float(line.strip().split('\t')[-1])
        if q_value >= 50:
            if q_value > gs_redundant_bin[gs_bin]:
                gs_redundant_bin[gs_bin]=q_value
                gs_redundant_bin2[gs_bin]=line
                try:
                    bin_with_microdiversity[gs_bin]+=1
                    i=bin_with_microdiversity[gs_bin]
                    bin_with_microdiversity2[gs_bin][i]=line
                except:
                    bin_with_microdiversity[gs_bin]=1
                    bin_with_microdiversity2[gs_bin]={}
                    bin_with_microdiversity2[gs_bin][1]=line

    f=open(test_bin_contigs+'_gs_bin_non-redundant.txt','w')
    f.write('Standard bin'+'\t'+'Test bin'+'\t'+'Standard bin length(bp)'+'\t'+'Test bin aligned length(bp)'+'\t'+'Standard bin aligned length(bp)'+'\t'+'Test bin aligned length(bp)'+'\t'+'Standard bin aligned ratio(%)'+'\t'+'Test bin aligned ratio(%)'+'\t'+'Contaminated ratio(%)'+'\t'+'CPN-CTN(%)'+'\t'+'CPN-3*CTN(%)'+'\t'+'CPN-5*CTN(%)'+'\n')
    for gs_bin in gs_redundant_bin2.keys():
        f.write(gs_redundant_bin2[gs_bin])
    f.close()

    test_bin_micro, xn={}, 0
    for line in open(str(test_bin_contigs)+'_bin_similarity_contamination_with_standard.txt','r'):
        xn+=1
        if xn >= 2:
            gs_bin=str(line).strip().split('\t')[0]
            test_bin=str(line).strip().split('\t')[1]
            q_value=float(str(line).strip().split('\t')[-1])
            if q_value >= 50:
                try:
                    test_bin_micro[test_bin][gs_bin]=line
                except:
                    test_bin_micro[test_bin]={}
                    test_bin_micro[test_bin][gs_bin]=line

    try:
        os.mkdir(test_bin_contigs+'_temp_orf_folder')
    except:
        print(test_bin_contigs+'_temp_orf_folder already existed. Recreated.')
        # os.system('rm -rf '+test_bin_contigs+'_temp_orf_folder')
        # os.mkdir(test_bin_contigs+'_temp_orf_folder')

    passed, num_orfs, orfs_len={}, {}, {}
    f=open(test_bin_contigs+'_with_potential_microdiversity.txt','w')
    f.write('Standard bin'+'\t'+'Test bin'+'\t'+'Standard bin length(bp)'+'\t'+'Test bin aligned length(bp)'+'\t'+'Standard bin aligned length(bp)'+'\t'+'Test bin aligned length(bp)'+'\t'+'Standard bin aligned ratio(%)'+'\t'+'Test bin aligned ratio(%)'+'\t'+'Contaminated ratio(%)'+'\t'+'CPN-CTN(%)'+'\t'+'CPN-3*CTN(%)'+'\t'+'CPN-5*CTN(%)'+'\n')
    for test_bin in test_bin_micro.keys():
        if len(test_bin_micro[test_bin]) >= 2:
            print(str(test_bin)+' '+str(test_bin_micro[test_bin]))
            os.chdir(pwd+'/'+test_bin_folder)



            os.system('prodigal -i '+str(test_bin)+' -d '+str(test_bin)+'.orfs.fna -m -p meta -q')
            orfs_len[str(test_bin)+'.orfs.fna']={}
            os.system('mv '+str(test_bin)+'.orfs.fna '+pwd+'/'+test_bin_contigs+'_temp_orf_folder')
            os.chdir(pwd)
            
            for gs_bin in test_bin_micro[test_bin].keys():
                f.write(str(test_bin_micro[test_bin][gs_bin]))
                if gs_bin not in passed.keys():
                    os.chdir(pwd+'/'+standard_bin_folder)
                    os.system('prodigal -i '+str(gs_bin)+' -d '+str(gs_bin)+'.orfs.fna -m -p meta -q')
                    orfs_len[str(gs_bin)+'.orfs.fna']={}
                    os.system('mv '+str(gs_bin)+'.orfs.fna '+pwd+'/'+test_bin_contigs+'_temp_orf_folder')
                    passed[gs_bin]=0
            os.chdir(pwd)
    f.close()

    os.chdir(pwd+'/'+test_bin_contigs+'_temp_orf_folder')
    for bins in orfs_len.keys():
        num=0
        for record in SeqIO.parse(bins,'fasta'):
            num+=1
            orfs_len[bins][record.id]=len(record.seq)
        num_orfs[bins]=num
    os.chdir(pwd)

    try:
        os.system('mkdir '+test_bin_contigs+'_summary')
    except:
        xxx=1
    os.system('mv *_basic_information_of_bins.txt '+str(test_bin_contigs)+'_bin_similarity_contamination_with_standard.txt '+test_bin_contigs+'_test_non-redundant_bin_similarity_contamination_with_standard.txt '+test_bin_contigs+'_gs_bin_non-redundant.txt '+test_bin_contigs+'_with_potential_microdiversity.txt '+test_bin_contigs+'_summary')

    aligned_orf, aligned_orf_gs_hit, aligned_orf_gs_hit2, aligned_orf2, aligned_orf3, aligned_orf_gs_bin, aligned_orf_gs_bin2, aligned_orf_gs_bin_x ={}, {}, {}, {}, {}, {}, {}, {}
    os.chdir(pwd+'/'+test_bin_contigs+'_temp_orf_folder')
    for test_bin in test_bin_micro.keys():
        aligned_orf_gs_bin[test_bin], aligned_orf_gs_bin_x[test_bin] ={}, {}
        if len(test_bin_micro[test_bin]) >= 2:
            os.system('makeblastdb -in '+str(test_bin)+'.orfs.fna -dbtype nucl')
            for gs_bin in test_bin_micro[test_bin].keys():
                os.system('blastn -db '+str(test_bin)+'.orfs.fna -query '+str(gs_bin)+'.orfs.fna -evalue 1e-20 -outfmt 6 -num_threads '+str(num_threads)+' -out '+str(test_bin)+'_vs_'+(gs_bin)+'.txt')
                fx=open('Filtated_'+str(test_bin)+'_vs_'+(gs_bin)+'.txt','w')
                aligned_orf[str(test_bin)+'\t'+str(gs_bin)]=0
                aligned_orf2[str(test_bin)+'\t'+str(gs_bin)]={}
                aligned_orf3[str(test_bin)+'\t'+str(gs_bin)]={}
                aligned_orf_gs_hit[str(test_bin)+'\t'+str(gs_bin)]=0
                aligned_orf_gs_hit2[str(test_bin)+'\t'+str(gs_bin)]={}
                for line in open(str(test_bin)+'_vs_'+(gs_bin)+'.txt','r'):
                    gs_bin_orf=str(line).strip().split('\t')[0].strip()
                    test_bin_orf=str(line).strip().split('\t')[1].strip()
                    sim=float(str(line).strip().split('\t')[2].strip())
                    ali_seq=int(str(line).strip().split('\t')[3].strip())
                    gs_bin_orf_s=int(str(line).strip().split('\t')[6].strip())
                    gs_bin_orf_e=int(str(line).strip().split('\t')[7].strip())
                    gs_ali_len=abs(gs_bin_orf_e-gs_bin_orf_s)
                    test_bin_orf_s=int(str(line).strip().split('\t')[8].strip())
                    test_bin_orf_e=int(str(line).strip().split('\t')[9].strip())
                    test_ali_len=abs(test_bin_orf_e-test_bin_orf_s)
                    gs_orf_ali_ratio=100*gs_ali_len/orfs_len[str(gs_bin)+'.orfs.fna'][gs_bin_orf]
                    test_orf_ali_ratio=100*test_ali_len/orfs_len[str(test_bin)+'.orfs.fna'][test_bin_orf]
                    if sim >= 99 and ali_seq >= 100:
                        try:
                            aligned_orf_gs_bin[test_bin][str(gs_bin)].append(test_bin_orf)
                            aligned_orf_gs_bin_x[test_bin][test_bin_orf][str(gs_bin)]=0
                        except:
                            aligned_orf_gs_bin[test_bin][str(gs_bin)]=[]
                            aligned_orf_gs_bin_x[test_bin][test_bin_orf]={}
                            aligned_orf_gs_bin[test_bin][str(gs_bin)].append(test_bin_orf)
                            aligned_orf_gs_bin_x[test_bin][test_bin_orf][str(gs_bin)]=0

                        fx.write(line)
                        try:
                            aligned_orf2[str(test_bin)+'\t'+(gs_bin)][gs_bin_orf]+=1
                        except:
                            aligned_orf2[str(test_bin)+'\t'+(gs_bin)][gs_bin_orf]=1
                        
                        try:
                            aligned_orf3[str(test_bin)+'\t'+str(gs_bin)][test_bin_orf]+=1
                        except:
                            aligned_orf[str(test_bin)+'\t'+(gs_bin)]+=1
                            aligned_orf3[str(test_bin)+'\t'+str(gs_bin)][test_bin_orf]=0
                        
                        try:
                            aligned_orf_gs_hit2[str(test_bin)+'\t'+str(gs_bin)][gs_bin_orf]+=1
                        except:
                            aligned_orf_gs_hit[str(test_bin)+'\t'+str(gs_bin)]+=1
                            aligned_orf_gs_hit2[str(test_bin)+'\t'+str(gs_bin)][gs_bin_orf]=0

                fx.close()
    # os.system('rm *.nin *.nhr *.nsq')

    print('Finding common and uniq ORFs')
    aligned_orf_gs_bin_x_level, hit_1_orf_source_bin ={}, {}
    for test_bin in aligned_orf_gs_bin_x.keys():
        aligned_orf_gs_bin_x_level[test_bin]={}
        hit_1_orf_source_bin[test_bin]={}
        for orfs in aligned_orf_gs_bin_x[test_bin].keys():
            num=len(aligned_orf_gs_bin_x[test_bin][orfs])
            try:
                aligned_orf_gs_bin_x_level[test_bin][num]+=1
            except:
                aligned_orf_gs_bin_x_level[test_bin][num]=1
            if num == 1:
                hit_1_orf_source_bin[test_bin][orfs]=str(aligned_orf_gs_bin_x[test_bin][orfs]).split(':')[0].replace('{','').replace('\'','').strip()
    
    hit_1_orf_source_bin2={}
    for test_bin in hit_1_orf_source_bin.keys():
        hit_1_orf_source_bin2[test_bin]={}
        for orfs in hit_1_orf_source_bin[test_bin].keys():
            gs_bin=hit_1_orf_source_bin[test_bin][orfs]
            try:
                hit_1_orf_source_bin2[test_bin][gs_bin]+=1
            except:
                hit_1_orf_source_bin2[test_bin][gs_bin]=1

    os.chdir(pwd+'/'+test_bin_contigs+'_summary')
    f=open(str(test_bin_contigs)+'_aligned_orfs_sumamry.txt','w')
    f.write('Test_bin'+'\t'+'NO. ORFs in test bin'+'\t'+'NO. ORFs hit in test bin'+'\t'+'Source of 1-hit only ORFs'+'\n')
    for test_bin in aligned_orf_gs_bin_x_level.keys():
        if len(aligned_orf_gs_bin[test_bin]) >= 2:
            test_bin_orf=test_bin+'.orfs.fna'
            f.write(str(test_bin)+'\t'+str(num_orfs[test_bin_orf])+'\t'+str(aligned_orf_gs_bin_x_level[test_bin]))
            if test_bin in hit_1_orf_source_bin2.keys():
                for gs_bin in hit_1_orf_source_bin2[test_bin].keys():
                    f.write('\t'+str(gs_bin)+': '+str(hit_1_orf_source_bin2[test_bin][gs_bin]))
            f.write('\n')
    f.close()

    f=open(str(test_bin_contigs)+'_aligned_orfs_sumamry2.txt','w')
    f.write('Test_bin'+'\t'+'GS_bin'+'\t'+'NO. aligned ORFs in test_bin'+'\t'+'NO. ORFs in test_bin'+'\t'+'Ratio aligned ORFs in test_bin'+'\t'+'NO. aligned ORFs in GS_bin'+'\t'+'NO. ORFs in GS_bin'+'\t'+'Ratio aligned ORFs in GS_bin'+'\t'+'No. Common_ORFs'+'\t'+'NO. Uniq ORFs in GS_bin'+'\n')
    for item in aligned_orf.keys():
        test_bin=item.split('\t')[0]+'.orfs.fna'
        gs_bin=item.split('\t')[1]+'.orfs.fna'
        orf_num=aligned_orf[item]
        test_bin_orfs_num=num_orfs[test_bin]
        ts_ratio=100*orf_num/test_bin_orfs_num
        gs_hit_orf_num=aligned_orf_gs_hit[item]
        gs_bin_orfs_num=num_orfs[gs_bin]
        gs_ratio=100*gs_hit_orf_num/gs_bin_orfs_num
        f.write(item+'\t'+str(orf_num)+'\t'+str(test_bin_orfs_num)+'\t'+str(ts_ratio)+'\t'+str(gs_hit_orf_num)+'\t'+str(gs_bin_orfs_num)+'\t'+str(gs_ratio)+'\n')
    f.close()
    # os.system('mv '+str(test_bin_contigs)+'_aligned_orfs_sumamry.txt '+pwd+'/'+test_bin_contigs+'_summary')
    os.chdir(pwd)
    # os.system('rm *.nin *.nhr *.nsq')
    print('Done')

if __name__ == '__main__': 
    pwd=os.getcwd()
    num_threads=60
    standard_bin_folder='CAMI_high_standard_genomes'
    test_bin_folder='BestBinset_outlier_refined_filtrated_retrieved_re-assembly_OLC_500-99-90'
    bin_length={}
    standard_bin=cat_contigs(standard_bin_folder, pwd)
    standard_bin_contigs=standard_bin[0]
    bin_length.update(standard_bin[1])
    test_bin=cat_contigs(test_bin_folder, pwd)
    test_bin_contigs=test_bin[0]
    bin_length.update(test_bin[1])

    Contigs_aligner(standard_bin_contigs, test_bin_contigs, bin_length, standard_bin_folder, test_bin_folder, num_threads, pwd)
