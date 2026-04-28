#!/usr/bin/env python

"""
Helper utilities for cleaning up intermediate BASALT files.

Implemented in pure Python so missing files / unmatched globs do not
spam stderr — all operations are best-effort and silent on no-op.
"""

import glob
import os
import shutil
import tarfile


def _rm(*patterns):
    """Remove files / directories matching any of the given glob patterns."""
    for pattern in patterns:
        for path in glob.glob(pattern):
            try:
                if os.path.isdir(path) and not os.path.islink(path):
                    shutil.rmtree(path, ignore_errors=True)
                else:
                    os.remove(path)
            except OSError:
                pass


def _move_into(dst, *patterns):
    """Move every match of *patterns* into directory *dst*. Returns True iff
    at least one file was moved (so the caller can decide whether to archive)."""
    moved = False
    for pattern in patterns:
        for path in glob.glob(pattern):
            if not moved:
                os.makedirs(dst, exist_ok=True)
            try:
                shutil.move(path, dst)
                moved = True
            except (OSError, shutil.Error):
                pass
    return moved


def _tar_gz(archive, *patterns):
    """Create *archive* (.tar.gz) containing everything matching *patterns*.
    No-op when nothing matches."""
    members = []
    for pattern in patterns:
        members.extend(glob.glob(pattern))
    if not members:
        return
    try:
        with tarfile.open(archive, 'w:gz') as tf:
            for member in members:
                tf.add(member)
    except (OSError, tarfile.TarError):
        pass


def cleanup(assembly_list):
    """
    Remove intermediate index/coverage files and archive key matrices.

    Parameters
    ----------
    assembly_list : list
        List of assembly file names; their numbered copies (``1_<name>``,
        ``2_<name>`` ...) are removed.
    """
    # BLAST index leftovers
    _rm('*.njs', '*.ndb', '*.nto', '*.ntf', '*.not', '*.nos')

    # Coverage / depth / connection matrices → archive
    backup = 'Coverage_depth_connection_SimilarBin_files_backup'
    if _move_into(backup,
                  '*.depth.txt', 'Coverage_matrix_*', 'Combat_*',
                  'condense_connections_*', 'Connections_*', 'Similar_bins.txt'):
        _tar_gz(backup + '.tar.gz', backup)
    if os.path.isdir(backup):
        shutil.rmtree(backup, ignore_errors=True)

    # Bulky scratch directories
    _rm('*_kmer', 'bin_coverage', 'Bin_coverage_after_contamination_removal',
        'bin_comparison_folder', 'bin_extract-eleminated-selected_contig',
        'Bins_blast_output')

    # Per-group archives
    _tar_gz('Group_comparison_files.tar.gz', '*_comparison_files')
    _tar_gz('Group_Bestbinset.tar.gz', '*_BestBinsSet')
    _tar_gz('Group_genomes.tar.gz', '*_genomes')
    _rm('*_sr_bins_seq')
    _tar_gz('Binsets_backup.tar.gz', 'BestBinse*')

    # Folders that have been archived above (or are otherwise discardable)
    _rm('*_comparison_files', '*_checkm', '*_genomes', '*_BestBinset',
        '*_BestBinsSet', 'BestBinse*', 'Deep_retrieved_bins',
        'coverage_deep_refined_bins', 'S6_coverage_filtration_matrix',
        'S6_TNF_filtration_matrix', 'split_blast_output',
        'TNFs_deep_refined_bins')

    _rm('*_checkpoint.txt')
    _rm('Merged_seqs_*')
    _rm('*.bt2', 'Outlier_in_threshold*', 'Summary_threshold*',
        'Refined_total_bins_contigs.fa', 'Total_bins.fa')

    for i in range(1, 20):
        _rm('*_deep_retrieval_' + str(i))

    _rm('*_MP_1', '*_MP_2', '*_gf_lr_polished', '*_gf_lr', '*_gf_lr_mod',
        '*_gf_lr_checkm', '*_long_read')

    # Misc per-run logs / scratch files
    _rm('Bin_reads_summary.txt', 'Depth_total.txt', 'Basalt_log.txt',
        'Assembly_mo_list.txt', 'Assembly_MoDict.txt', '*_gf_lr_blasted.txt',
        'Bestbinset_list.txt', 'Bin_extract_contigs_after_coverage_filtration.txt',
        'Bin_lw.txt', 'Bin_record_error.txt')
    _rm('Bins_folder.txt', 'BLAST_output_error.txt', 'Concoct_*',
        'condensed.cytoscape*', 'cytoscape.*')
    _rm('Hybrid_re-assembly_status.txt', 'Mapping_log_*',
        'OLC_merged_error_blast_results.txt',
        'Potential_contaminted_seq_vari.txt',
        'Reassembled_bins_comparison.txt', 'Rejudge_clean.txt')
    _rm('Remained_seq*', 'Remapped_depth_test.txt', 'Re-mapped_depth.txt',
        'Remapping.fasta', 'TNFs_exceptional_contigs.txt',
        'Total_contigs_after_OLC_reassembly.fa')
    _rm('PE_r1_*', 'PE_r2_*')

    # Numbered assembly copies (1_<name>, 2_<name>, ...)
    for i, name in enumerate(assembly_list, start=1):
        _rm('{}_{}'.format(i, name))


if __name__ == '__main__':
    assembly_list = ['8_medium_S001_SPAdes_scaffolds.fasta',
                     '10_medium_cat_SPAdes_scaffolds.fasta']
    cleanup(assembly_list)
