#!/usr/bin/env python

"""
Step S9: Short-read-based reassembly of bins.

Unified entry for CheckM2 and CheckM backends. Call ``set_qc_backend(name)``
before invoking ``re_assembly_main`` to select the QC tool.
"""

from Bio import SeqIO
import sys, os, threading, copy, math
import concurrent.futures
from multiprocessing import Pool

from basalt.qc_backend import get_backend


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
            'S9_Reassembly: QC backend not set. '
            "Call set_qc_backend('checkm2' or 'checkm') first."
        )
    return _BACKEND


def mod_bin(binset_folder):
    """
    Create a re-indexed version of a binset with sequential bin IDs.

    Parameters
    ----------
    binset_folder : str
        Folder containing the original per-bin FASTA and a QC results
        TSV produced by either CheckM2 (``quality_report.tsv``) or
        CheckM (``storage/bin_stats_ext.tsv``).

    Returns
    -------
    str
        Name of the new binset folder (``<binset_folder>_mod``).
    str
        Path to the combined ``Total_bins.fa`` in the working directory.
    list of str
        List of new bin IDs (``bin1``, ``bin2``, ...).
    dict
        Mapping ``bin_id -> quality_metrics`` (unified schema, see
        ``qc_backend`` module docstring).
    """
    backend = _require_backend()
    pwd = os.getcwd()
    bins_checkm = {}
    try:
        os.mkdir(str(binset_folder)+'_mod')
    except:
        print(str(binset_folder)+'_mod is exist. Re-create the folder')
        os.system('rm -rf '+str(binset_folder)+'_mod')
        os.mkdir(str(binset_folder)+'_mod')

    os.chdir(pwd+'/'+binset_folder)
    f=open('Mod_contig_id.txt','w')
    f1=open('Bin_name_mod.txt','w')
    mod_tsv_name = 'Bin_name_mod_'+backend.results_filename
    f2=open(mod_tsv_name,'w')
    f2.write('Bin_ID'+'\t'+'Genome_size'+'\t'+'Completeness'+'\t'+'Contamination'+'\t'+backend.contig_size_key+'\n')
    n, record_seq, mod_bin_list, mod_bin_dict, bin_contig_num = 0, {}, [], {}, {}
    for root, dirs, files in os.walk(pwd+'/'+binset_folder):
        for file in files:
            hz=str(file).split('.')[-1]
            if 'fa' in hz:
                n+=1
                m=0
                record_seq['bin'+str(n)]={}
                mod_bin_list.append('bin'+str(n))
                checkm_name_list=str(file).split('.')
                checkm_name_list.remove(checkm_name_list[-1])
                checkm_name='.'.join(checkm_name_list)
                mod_bin_dict[checkm_name]='bin'+str(n)
                f1.write(str(file)+'\t'+'bin'+str(n)+'.fa'+'\n')
                for record in SeqIO.parse(file, 'fasta'):
                    m+=1
                    f.write('bin'+str(n)+'_'+str(m)+'\t'+str(record.id)+'\n')
                    record_seq['bin'+str(n)]['bin'+str(n)+'_'+str(m)]=str(record.seq)
                bin_contig_num['bin'+str(n)]=m

    if n == 0:
        f.close()
        f1.close()
        f2.close()
        os.chdir(pwd)
        raise RuntimeError(
            'S9_Reassembly: binset {!r} selected for reassembly contains no '
            '.fa bins. An upstream step (outlier removal / filtration / OLC) '
            'produced an empty binset, so Total_bins.fa would be empty. '
            'Aborting before the Bowtie2/mapping stage. Check the folder named '
            'on the "8th contig OLC" line of Basalt_checkpoint.txt and the '
            'outlier/filtration outputs upstream of it.'.format(binset_folder)
        )

    raw_metrics = backend.parse_results(pwd+'/'+binset_folder)
    for original_name, metrics in raw_metrics.items():
        if original_name not in mod_bin_dict:
            continue
        mod_name = mod_bin_dict[original_name]
        bins_checkm[mod_name] = {
            'Completeness': float(metrics.get('Completeness', 0.0)),
            'Contamination': float(metrics.get('Contamination', 0.0)),
            'Genome size': int(metrics.get('Genome size', 0)),
            'contig_size': float(metrics.get('contig_size', 0.0)),
            'contig_size_key': metrics.get('contig_size_key', backend.contig_size_key),
        }
        f2.write(
            mod_name+'\t'
            +str(bins_checkm[mod_name]['Genome size'])+'\t'
            +str(bins_checkm[mod_name]['Completeness'])+'\t'
            +str(bins_checkm[mod_name]['Contamination'])+'\t'
            +str(bins_checkm[mod_name]['contig_size'])+'\n'
        )
    f.close()
    f1.close()
    f2.close()

    os.system('mv Mod_contig_id.txt Bin_name_mod.txt '+mod_tsv_name+' '+pwd+'/'+str(binset_folder)+'_mod')

    os.chdir(pwd+'/'+str(binset_folder)+'_mod')
    f1=open('Total_bins.fa','w')
    for bin_id in record_seq.keys():
        f=open(str(bin_id)+'.fa','w')
        for contigs in record_seq[bin_id].keys():
            f.write('>'+str(contigs)+'\n'+str(record_seq[bin_id][contigs])+'\n')
            f1.write('>'+str(contigs)+'\n'+str(record_seq[bin_id][contigs])+'\n')
        f.close()
    f1.close()
    os.system('mv Total_bins.fa '+pwd)
    os.chdir(pwd)
    return str(binset_folder)+'_mod', 'Total_bins.fa', mod_bin_list, bins_checkm


def parse_sam(sam_file, fq, pair, n, num_threads=4):
    """
    Parse a Bowtie2 SAM file and split paired alignments by bin.

    Single-pass implementation: caches per-(bin, read_id_name) alignments
    until a second hit for the same pair arrives, then queues both records
    into a per-output-file write buffer. The buffers are flushed once at
    the end with one open/close per FASTQ file, parallelised across files
    with a thread pool. Updates ``fq`` and ``pair`` in place to preserve
    the original API contract.
    """
    print('Parsing', sam_file)
    f_not_mapped_reads = open('Not_mapped_reads.txt', 'a')

    # bin_id -> {rid: [(read_num, record), ...]}  alignments awaiting a 2nd hit
    pending = {}
    # (bin_id, read_num) -> [record, ...]   write buffers, one entry per output fastq
    buffers = {}

    m = 0
    for line in open(sam_file, 'r'):
        m += 1
        flist = line.split('\t')
        if len(flist) < 12:
            if m % 1000000 == 0:
                print('Parsed', m, 'lines')
            continue
        bin_id = flist[2].split('_')[0]
        read_id = flist[0]
        rid = str(n) + '_' + read_id.split('_')[0]
        read_num = (str(n) + '_' + read_id).split('_')[-1]

        b_pending = pending.setdefault(bin_id, {})
        b_fq = fq.setdefault(bin_id, {})

        if rid in b_fq:
            # Already emitted both mates for this rid → multimapper, drop & log
            f_not_mapped_reads.write(rid + '\n')
        else:
            fq_seq = flist[9] + '\n+\n' + flist[10] + '\n'
            rec = '@' + (str(n) + '_' + read_id)[:-2] + ' ' + read_num + '\n' + fq_seq
            b_pending.setdefault(rid, []).append((read_num, rec))
            if len(b_pending[rid]) == 2:
                for rn, r in b_pending[rid]:
                    buffers.setdefault((bin_id, rn), []).append(r)
                b_fq[rid] = 2
                del b_pending[rid]

        if m % 1000000 == 0:
            print('Parsed', m, 'lines')

    # rids that never reached 2 hits → mirror legacy ``pair`` state and log
    for bin_id, b_pending in pending.items():
        b_pair = pair.setdefault(bin_id, {})
        for rid in b_pending:
            b_pair[rid] = 1
            f_not_mapped_reads.write(rid + '\n')
    f_not_mapped_reads.close()

    f_summary = open('Bin_reads_summary.txt', 'w')
    for bin_id in fq.keys():
        f_summary.write(str(bin_id) + ' SEQ number:' + str(len(fq[bin_id])) + '\n')
    f_summary.close()

    if buffers:
        print('Flushing FASTQ buffers (' + str(len(buffers)) + ' files, '
              + str(max(1, int(num_threads))) + ' threads)')

        def _flush(item):
            (bin_id, rn), recs = item
            with open(str(bin_id) + '_seq_R' + str(rn) + '.fq', 'a') as f:
                f.writelines(recs)

        with concurrent.futures.ThreadPoolExecutor(
                max_workers=max(1, int(num_threads))) as pool:
            for _ in pool.map(_flush, list(buffers.items())):
                pass


def parse_lr_sam(sam_file, long_read, sn):
    """Group long reads by bin based on SAM alignments."""
    print('Reading long reads id '+str(long_read))
    bin_lr, bin_lr2, lr_bin, lr_bin2 = {}, {}, {}, {}
    m, m1, m2 = 0, 0, 0
    for line in open(sam_file,'r'):
        m1+=1
        flist=str(line).split('\t')
        if len(flist) >= 12:
            read_id=flist[0]
            bin_id=flist[2].split('_')[0]
            bin_lr[bin_id]={}
            lr_bin[read_id]={}

        if m1 % 1000000 == 0:
            print('Read', m1,'lines')

    print('Collecting long reads '+str(long_read))
    for line in open(sam_file,'r'):
        m2+=1
        flist=str(line).split('\t')
        if len(flist) >= 12:
            bin_id=flist[2].split('_')[0]
            seqs=flist[9]
            if 'bin' in bin_id and len(seqs) >= 100:
                read_id=flist[0]
                bin_lr[bin_id][read_id]=''
                lr_bin[read_id][bin_id]=''

        if m2 % 1000000 == 0:
            print('Read', m2,'lines')

    bin_lr2=copy.deepcopy(bin_lr)
    lr_bin2=copy.deepcopy(lr_bin)

    for item in bin_lr2.keys():
        if len(bin_lr2[item]) == 0:
            del bin_lr[item]

    for item in lr_bin2.keys():
        if len(lr_bin2[item]) == 0:
            del lr_bin[item]

    f=open('Bin_long_read'+str(sn)+'.txt','w')
    for bin_id in bin_lr.keys():
        f.write(str(bin_id)+'\t'+str(bin_lr[bin_id])+'\n')
        try:
            f1=open(str(bin_id)+'_lr.fq','a')
        except:
            f1=open(str(bin_id)+'_lr.fq','w')
        f1.close()
    f.close()

    f=open('Long_read_bin'+str(sn)+'.txt','w')
    for lr in lr_bin.keys():
        f.write(str(lr)+'\t'+str(lr_bin[lr])+'\n')
    f.close()

    n, record_bin_line, record_bin_line2 = 0, {}, {}
    for line in open(long_read,'r'):
        n+=1
        record_bin_line[n]=[]

    print('Splitting long reads to different bins '+str(long_read))
    n = 0
    for line in open(long_read,'r'):
        n+=1
        m=n-1
        if m % 4 == 0:
            seq_id=str(line).strip().split(' ')[0].split('@')[1]
            if seq_id in lr_bin.keys():
                for bin_id in lr_bin[seq_id].keys():
                    record_bin_line[n].append(bin_id)
                    record_bin_line[n+1].append(bin_id)
                    record_bin_line[n+2].append(bin_id)
                    record_bin_line[n+3].append(bin_id)

    record_bin_line2=copy.deepcopy(record_bin_line)
    for line in record_bin_line2.keys():
        if len(record_bin_line2[line]) == 0:
            del record_bin_line[line]
    seq_num=int(len(record_bin_line)/4)
    print(str(seq_num)+' reads from '+str(long_read)+' will be splitted into different bins')

    n1=0
    for line in open(long_read,'r'):
        n1+=1
        if n1 in record_bin_line.keys():
            for bin_id in record_bin_line[n1]:
                f1=open(str(bin_id)+'_lr.fq','a')
                f1.write(line)
                f1.close()
    if n1 % 1000000 == 0:
        print('Parse', n1,'lines')
    print('Long reads '+str(long_read)+' splitting done!')
    return bin_lr


def mapping_sr(total_fa, datasets_list, fq, pair, num_threads):
    """Map short reads to bins and prepare bin-specific FASTQ files."""
    os.system('bowtie2-build '+str(total_fa)+' '+str(total_fa))
    n = 0
    for item in datasets_list.keys():
        n+=1
        os.system('bowtie2 -p '+str(num_threads)+' -x '+str(total_fa)+' -1 '+str(datasets_list[item][0])+' -2 '+str(datasets_list[item][1])+' -S '+str(item)+'.sam -q --no-unal')
        parse_sam(str(item)+'.sam', fq, pair, n, num_threads=num_threads)
        os.system('rm '+str(item)+'.sam')


def parse_checkm(checkm_containing_folder, pwd):
    """
    Parse QC reports produced by either CheckM2 or CheckM for S9.

    Always call with the top-level output folder. The backend walks the
    tree, so CheckM's ``storage/bin_stats_ext.tsv`` is found without the
    caller having to append ``/storage``.
    """
    backend = _require_backend()
    return backend.parse_results(os.path.join(pwd, checkm_containing_folder))


def reassembly(bin_seq, reassembly_bin_folder, num_threads, bins_seq_folder,
               long_read, ram, pwd, sensitive='sensitive'):
    """Per-bin short-read reassembly. Always runs MEGAHIT (--presets meta-sensitive); also runs SPAdes when sensitive=='more-sensitive'."""
    for item in bin_seq.keys():
        megahit_out=str(item)+'_megahit_reassembly'
        os.system('rm -rf '+megahit_out)
        os.system('megahit --presets meta-sensitive -1 '+str(bin_seq[item][0])+' -2 '+str(bin_seq[item][1])+' -o '+megahit_out+' --min-contig-len 1000 -t '+str(num_threads))
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
            os.system('spades.py -1 '+str(bin_seq[item][0])+' -2 '+str(bin_seq[item][1])+' -o '+str(item)+'_spades_reassembly --careful -t '+str(num_threads)+' -m '+str(ram))
            os.chdir(str(pwd)+'/'+str(item)+'_spades_reassembly')
            xxxx=0
            if os.path.isfile("contigs.fasta") == True:
                f=open(str(item)+'_SPAdes_re-assembly_contigs.fa','w')
                for record in SeqIO.parse('contigs.fasta', 'fasta'):
                    if len(record.seq) >= 1000:
                        xxxx+=1
                        f.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
                f.close()
            else:
                os.chdir(str(pwd))
                os.system('spades.py -1 '+str(bin_seq[item][0])+' -2 '+str(bin_seq[item][1])+' -o '+str(item)+'_spades_reassembly -t '+str(num_threads)+' -m '+str(ram))
                os.chdir(str(pwd)+'/'+str(item)+'_spades_reassembly')
                if os.path.isfile("contigs.fasta") == True:
                    f=open(str(item)+'_SPAdes_re-assembly_contigs.fa','w')
                    for record in SeqIO.parse('contigs.fasta', 'fasta'):
                        if len(record.seq) >= 1000:
                            xxxx+=1
                            f.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
                    f.close()
            os.chdir(str(pwd))

            if xxxx >= 2:
                os.system('mv '+str(item)+'_spades_reassembly/'+str(item)+'_SPAdes_re-assembly_contigs.fa '+str(reassembly_bin_folder)+'/'+str(item)+'_SPAdes_re-assembly_contigs.fa')
            os.system('rm -rf '+str(item)+'_spades_reassembly')

        os.system('mv '+str(bin_seq[item][0])+' '+str(bin_seq[item][1])+' '+str(bins_seq_folder))


def unicycler_mul(item, pwd, sr_folder, bin_seq, bin_lr, t_p_p, reassembly_bin_folder, bins_seq_folder):
    os.system('unicycler -1 '+str(pwd)+'/'+str(sr_folder)+'/'+str(bin_seq[item][0])+' -2 '+str(pwd)+'/'+str(sr_folder)+'/'+str(bin_seq[item][1])+' -l '+str(bin_lr[item])+' -o '+str(item)+'_unicycler_reassembly -t '+str(t_p_p)+' --mode conservative --no_pilon')
    os.system('mv '+str(item)+'_unicycler_reassembly/assembly.fasta '+str(reassembly_bin_folder)+'/'+str(item)+'_UNICYCLER_re-assembly_contigs.fa')
    os.system('mv '+str(bin_lr[item])+' '+str(pwd)+'/'+str(bins_seq_folder))
    os.system('rm -rf '+str(item)+'_unicycler_reassembly')


def reassembly_lr(bin_seq, bin_lr, reassembly_bin_folder, num_threads,
                  bins_seq_folder, sr_folder, ram, pwd):
    """Perform hybrid reassembly of bins using both short and long reads."""
    try:
        os.system('java -Xmx'+str(ram)+'G -jar pilon-1.23.jar')
    except:
        print('pilon running in default ram')

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

    if num_project == 1:
        for item in bin_lr.keys():
            if item in bin_seq.keys():
                print('Reassembling '+str(item))
                os.system('unicycler -1 '+str(pwd)+'/'+str(sr_folder)+'/'+str(bin_seq[item][0])+' -2 '+str(pwd)+'/'+str(sr_folder)+'/'+str(bin_seq[item][1])+' -l '+str(bin_lr[item])+' -o '+str(item)+'_unicycler_reassembly -t '+str(num_threads)+' --mode conservative --no_pilon')
                os.system('mv '+str(item)+'_unicycler_reassembly/assembly.fasta '+str(reassembly_bin_folder)+'/'+str(item)+'_UNICYCLER_re-assembly_contigs.fa')
                os.system('mv '+str(bin_lr[item])+' '+str(pwd)+'/'+str(bins_seq_folder))
                os.system('rm -rf '+str(item)+'_unicycler_reassembly')
            else:
                xyz=0
    else:
        t_p_p=math.ceil(num_threads/num_project)
        print('Hybrid reassembly parallelism: '+str(num_project)+' bins assembled in parallel x '+str(t_p_p)+' threads/bin '
              '(--threads='+str(num_threads)+', --ram='+str(ram)+'G; Unicycler will report '+str(t_p_p)+' threads each by design)')
        pool=Pool(processes=num_project)
        for item in bin_lr.keys():
            if item in bin_seq.keys():
                print('Reassembling '+str(item))
                pool.apply_async(unicycler_mul, args=(item, pwd, sr_folder, bin_seq, bin_lr, t_p_p, reassembly_bin_folder, bins_seq_folder))
        pool.close()
        pool.join()

        os.system('mv '+item+' '+str(bins_seq_folder))
    os.system('rm -rf '+str(item)+'_unicycler_reassembly')


def bin_comparison(paired_bins, bin_checkm):
    pwd=os.getcwd()
    f=open('Reassembled_bins_comparison.txt','w')
    best_bin, best_bin_checkm={}, {}
    for item in paired_bins.keys():
        best_bin_checkm_name_list=item.split('.')
        best_bin_checkm_name_list.remove(best_bin_checkm_name_list[-1])
        best_bin_checkm_name='.'.join(best_bin_checkm_name_list)
        f.write(str(item)+'\t'+str(bin_checkm[best_bin_checkm_name])+'\n')
        for item2 in paired_bins[item]:
            reass_bin_checkm_name_list=item2.split('.')
            reass_bin_checkm_name_list.remove(reass_bin_checkm_name_list[-1])
            reass_bin_checkm_name='.'.join(reass_bin_checkm_name_list)
            f.write(str(item2)+'\t'+str(bin_checkm[reass_bin_checkm_name])+'\n')
            best_bin_cpn=bin_checkm[best_bin_checkm_name]['Completeness']
            best_bin_ctn=bin_checkm[best_bin_checkm_name]['Contamination']
            best_bin_ml=bin_checkm[best_bin_checkm_name]['contig_size']
            reass_bin_cpn=bin_checkm[reass_bin_checkm_name]['Completeness']
            reass_bin_ctn=bin_checkm[reass_bin_checkm_name]['Contamination']
            reass_bin_ml=bin_checkm[reass_bin_checkm_name]['contig_size']

            delta_cpn_ctn_bestbin=float(best_bin_cpn)-5*float(best_bin_ctn)
            delta_cpn_ctn_reass_bin=float(reass_bin_cpn)-float(5*reass_bin_ctn)

            if '_SPAdes_' in reass_bin_checkm_name or '_UNICYCLER_' in reass_bin_checkm_name:
                if delta_cpn_ctn_bestbin > delta_cpn_ctn_reass_bin:
                    best_bin_checkm_name=best_bin_checkm_name
                elif delta_cpn_ctn_bestbin < delta_cpn_ctn_reass_bin:
                    best_bin_checkm_name=reass_bin_checkm_name
                elif delta_cpn_ctn_bestbin == delta_cpn_ctn_reass_bin:
                    if reass_bin_ml > best_bin_ml:
                        best_bin_checkm_name=reass_bin_checkm_name
                    else:
                        best_bin_checkm_name=best_bin_checkm_name
                else:
                    continue
            elif '_IDBA_' in reass_bin_checkm_name:
                delta_idba = delta_cpn_ctn_reass_bin-delta_cpn_ctn_bestbin
                if delta_idba < 3:
                    best_bin_checkm_name=best_bin_checkm_name
                elif delta_idba >= 3:
                    best_bin_checkm_name=reass_bin_checkm_name
                else:
                    continue

        best_bin[best_bin_checkm_name+'.fa']=best_bin_checkm_name
        best_bin_checkm[best_bin_checkm_name]=bin_checkm[best_bin_checkm_name].copy()
    f.close()
    return best_bin, best_bin_checkm


def re_assembly_main(binset_folder, datasets_list, long_read,
                     hybri_reassembly, ram, num_threads, sensitive='sensitive'):
    """
    Entry point for S9 short-read reassembly.

    Parameters
    ----------
    binset_folder : str
        Folder containing the starting binset.
    datasets_list : dict
        Mapping ``sample_id -> [r1_fastq, r2_fastq]`` with paired-end reads.
    long_read : list of str
        List of long-read FASTQ files used for optional hybrid reassembly.
    hybri_reassembly : {'y', 'n'}
        Flag controlling whether hybrid reassembly should be performed.
    ram : int
        Maximum RAM (in GB) available to assemblers.
    num_threads : int
        Number of threads to use for mapping and assembly.
    """
    backend = _require_backend()
    pwd=os.getcwd()
    try:
        f=open('Assembly_status.txt','a')
    except:
        f=open('Assembly_status.txt','w')
    f.close()

    try:
        os.mkdir(binset_folder+'_re-assembly_binset')
    except:
        print(binset_folder+'_re-assembly_binset exists')

    assembled_bins={}
    try:
        os.mkdir(binset_folder+'_sr_bins_seq')
    except:
        print(binset_folder+'_sr_bins_seq exists')
        for root, dirs, files in os.walk(pwd+'/'+binset_folder+'_sr_bins_seq'):
            for file in files:
                if '.fq' in file and '_seq_' in file and 'bin' in file:
                    assembled_bins[str(file).split('_seq_')[0].strip()]=0

    A=mod_bin(binset_folder)
    mod_bin_folder=A[0]
    total_fa=A[1]
    mod_bin_list=A[2]
    original_bins_checkm=A[3]
    fq, pair, bin_seq, bin_checkm={}, {}, {}, {}
    bin_checkm=original_bins_checkm.copy()
    for bin_id in mod_bin_list:
        fq[str(bin_id)]={}
        pair[str(bin_id)]={}
        bin_seq[str(bin_id)]=[]
        bin_seq[str(bin_id)].append(str(bin_id)+'_seq_R1.fq')
        bin_seq[str(bin_id)].append(str(bin_id)+'_seq_R2.fq')
        if len(assembled_bins) == 0:
            f1=open(str(bin_id)+'_seq_R1.fq','w')
            f2=open(str(bin_id)+'_seq_R2.fq','w')
            f1.close()
            f2.close()

    x=0
    for line in open('Assembly_status.txt','r'):
        if 'Short-read mapping done!' in line:
            x=1

    if x == 0:
        f_not_mapped_reads=open('Not_mapped_reads.txt','w')
        f_not_mapped_reads.close()

        mapping_sr(total_fa, datasets_list, fq, pair, num_threads)
        f=open('Assembly_status.txt','a')
        f.write('Short-read mapping done!'+'\n')
        f.close()

    bin_seq2={}
    for item in bin_seq:
        if item not in assembled_bins.keys():
            bin_seq2[item]=bin_seq[item]

    x=0
    for line in open('Assembly_status.txt','r'):
        if 'Short-read assembly done!' in line:
            x=1

    if x == 0:
        reassembly(bin_seq2, binset_folder+'_re-assembly_binset', num_threads, binset_folder+'_sr_bins_seq', long_read, ram, pwd, sensitive=sensitive)
        f=open('Assembly_status.txt','a')
        f.write('Short-read assembly done!'+'\n')
        f.close()

    if len(long_read) != 0 and hybri_reassembly == 'y':
        print('Reassemblying long reads')
        assembled_bins={}
        try:
            os.mkdir(binset_folder+'_lr_bins_seq')
        except:
            print(binset_folder+'_lr_bins_seq exists')
            for root, dirs, files in os.walk(pwd+'/'+binset_folder+'_lr_bins_seq'):
                for file in files:
                    if '_lr.fq' in file and 'bin' in file:
                        assembled_bins[str(file).split('_lr')[0].strip()]=0
        os.chdir(pwd)

        if len(assembled_bins) == 0:
            print('Mapping long reads')
            x=0
            for line in open('Assembly_status.txt','r'):
                if 'Long-read mapping done!' in line:
                    x=1

            if x == 0:
                n, bin_lr=0, {}
                for lrs in long_read:
                    n+=1
                    os.system('minimap2 -t '+str(num_threads)+' -ax map-ont Total_bins.fa '+str(lrs)+' > lr'+str(n)+'.sam')
                    print('Splitting long reads '+str(n))
                    bin_lr.update(parse_lr_sam('lr'+str(n)+'.sam', lrs, n))

                f=open('Assembly_status.txt','a')
                f.write('Long-read mapping done!'+'\n')
                f.close()

                for i in range(1, n+1):
                    os.system('rm lr'+str(i)+'.sam')

        x=0
        for line in open('Assembly_status.txt','r'):
            if 'Long-read assembly done!' in line:
                x=1

        if x == 0:
            n, bin_lr=0, {}
            for lrs in long_read:
                n+=1
                for line in open('Bin_long_read'+str(n)+'.txt','r'):
                    bin_id=str(line).strip().split('\t')[0]
                    bin_lr[bin_id]=0

            bin_lr2={}
            for item in bin_lr.keys():
                if item not in assembled_bins.keys():
                    bin_lr2[item]=str(item)+'_lr.fq'

            reassembly_lr(bin_seq, bin_lr2, binset_folder+'_re-assembly_binset', num_threads, binset_folder+'_lr_bins_seq', binset_folder+'_sr_bins_seq', ram, pwd)
            f=open('Assembly_status.txt','a')
            f.write('Long-read assembly done!'+'\n')
            f.close()

    x=0
    for line in open('Assembly_status.txt','r'):
        if 'Assembly quality evaluation done!' in line:
            x=1

    if x == 0:
        print('Checking reassembled bins')
        os.chdir(str(binset_folder)+'_re-assembly_binset')
        for root, dirs, files in os.walk(pwd+'/'+str(binset_folder)+'_re-assembly_binset'):
            for file in files:
                try:
                    ff=open(file+'_filtrated.fa','w')
                    for record in SeqIO.parse(file, 'fasta'):
                        if len(record.seq) != 0:
                            ff.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
                    ff.close()
                    os.system('mv '+str(file)+'_filtrated.fa '+str(file))
                except:
                    xyzzzz=0
        os.chdir(pwd)

        # CheckM has historically been flaky on the reassembled binset; if the
        # results TSV is missing after the first run we wipe the output dir
        # and retry once before giving up. CheckM2 is reliable but the same
        # guard is harmless.
        out_dir = str(binset_folder)+'_re-assembly_binset_checkm'
        in_dir = str(binset_folder)+'_re-assembly_binset'
        ok, attempt = False, 0
        while not ok and attempt < 2:
            attempt += 1
            if attempt > 1:
                os.system('rm -rf '+out_dir)
            os.system(backend.build_cmd(num_threads, in_dir, out_dir, ext='fa'))
            results_file = backend.results_path(pwd+'/'+out_dir)
            try:
                with open(results_file, 'r') as _rf:
                    for _line in _rf:
                        ok = True
                        break
            except (FileNotFoundError, OSError):
                ok = False

        f=open('Assembly_status.txt','a')
        if ok:
            f.write('Assembly quality evaluation done!'+'\n')
        else:
            f.write('Reassembly QC error! Abandon!'+'\n')
            print('Reassembly QC error! Abandon!')
        f.close()

    os.chdir(pwd)
    reassembly_bin_checkm=parse_checkm(binset_folder+'_re-assembly_binset_checkm', pwd)
    bin_checkm.update(reassembly_bin_checkm)

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
            elif '_UNICYCLER_re-assembly_contigs.fa' in file:
                n_re=0
                for line in open(file,'r'):
                    n_re+=1
                    if n_re == 2:
                        break
                if n_re == 2:
                    o_bin=str(file).split('_UNICYCLER_re-assembly_contigs.fa')[0]+'.fa'
                    if o_bin not in paired_bins.keys():
                        paired_bins[o_bin]=[]
                        paired_bins[o_bin].append(str(file))
                    else:
                        paired_bins[o_bin].append(str(file))
            elif '_FLYe_re-assembly_contigs.fa' in file:
                n_re=0
                for line in open(file,'r'):
                    n_re+=1
                    if n_re == 2:
                        break
                if n_re == 2:
                    o_bin=str(file).split('_FLYe_re-assembly_contigs.fa')[0]+'.fa'
                    if o_bin not in paired_bins.keys():
                        paired_bins[o_bin]=[]
                        paired_bins[o_bin].append(str(file))
                    else:
                        paired_bins[o_bin].append(str(file))
    os.chdir(pwd)

    best_bin_after_c=bin_comparison(paired_bins, bin_checkm)
    best_bin=best_bin_after_c[0].copy()
    best_bin_checkm=best_bin_after_c[1].copy()

    try:
        os.mkdir(binset_folder+'_re-assembly')
    except:
        print(binset_folder+'_re-assembly exists')

    selected_bins={}
    for item in best_bin.keys():
        if '_re-assembly_contigs.fa' in item:
            s_name=item.split('_')[0]+'.fa'
            selected_bins[s_name]=1
            os.system('cp '+pwd+'/'+binset_folder+'_re-assembly_binset/'+item+' '+pwd+'/'+binset_folder+'_re-assembly')
        else:
            selected_bins[item]=1
            os.system('cp '+pwd+'/'+str(binset_folder)+'_mod/'+item+' '+pwd+'/'+binset_folder+'_re-assembly')

    os.chdir(pwd+'/'+str(binset_folder)+'_mod')
    for root, dirs, files in os.walk(pwd+'/'+str(binset_folder)+'_mod'):
        for file in files:
            if '.fa' in file:
                if file not in selected_bins.keys():
                    os.system('cp '+pwd+'/'+str(binset_folder)+'_mod/'+file+' '+pwd+'/'+binset_folder+'_re-assembly')
                    item_checkm_name=file.split('.fa')[0]
                    best_bin_checkm[item_checkm_name]=bin_checkm[item_checkm_name]
    os.chdir(pwd)

    os.chdir(pwd+'/'+binset_folder+'_re-assembly')
    f=open('Best_binset_after_re-assembly_quality_report.tsv','w')
    f.write('Bin_ID'+'\t'+'Genome_size'+'\t'+'Completeness'+'\t'+'Contamination'+'\t'+backend.contig_size_key+'\n')
    for item in best_bin_checkm.keys():
        f.write(
            str(item)+'\t'
            +str(best_bin_checkm[item].get('Genome size', 0))+'\t'
            +str(best_bin_checkm[item].get('Completeness', 0.0))+'\t'
            +str(best_bin_checkm[item].get('Contamination', 0.0))+'\t'
            +str(best_bin_checkm[item].get('contig_size', 0.0))+'\n'
        )
    f.close()
    os.chdir(pwd)

    print('Re-assembly done!')


if __name__ == '__main__':
    binset_folder='1_Opera_unpolished_cat_contigs.fasta_BestBinsSet_outlier_refined_filtrated_retrieved'
    datasets_list={'1':['PE_r1_RH_S001_insert_270_mate1.fq','PE_r2_RH_S001_insert_270_mate2.fq']}
    long_read=[]
    hybri_reassembly='n'
    num_threads=20
    ram=250
    set_qc_backend('checkm2')
    re_assembly_main(binset_folder, datasets_list, long_read, hybri_reassembly, ram, num_threads)
