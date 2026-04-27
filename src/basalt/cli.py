#!/usr/bin/python
# -*- coding: UTF-8 -*-

"""
Command-line interface for the BASALT metagenomic binning pipeline.

This script parses user arguments and dispatches to the corresponding
pipeline modules (autobinning, refinement, reassembly, data feeding,
and dereplication) depending on the selected mode and options.
"""

import time
import sys
import os
import argparse
import warnings
from glob import glob

try:
    from sklearn.exceptions import EfficiencyWarning
    warnings.filterwarnings("ignore", category=EfficiencyWarning)
except ImportError:
    pass


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='BASALT')
parser.add_argument('-a', '--assemblies', type=str, dest='assemblies',
                    help='List of assemblies, e.g.: as1.fa,as2.fa')
parser.add_argument('-s', '--shortreads', type=str, dest='sr_datasets',
                    help='List of paired-end reads, e.g.: r1_1.fq,r1_2.fq/r2_1.fq,r2_2.fq (paried_ends reads need \'/\' to seperate)')
parser.add_argument('-l', '--longreads', type=str, dest='long_datasets',
                    help='Including ont and pb dataset, excluding hifi dataset. List of long reads, e.g.: lr1.fq,lr2.fq')
parser.add_argument('-hf', '--hifi', type=str, dest='hifi_datasets',
                    help='Hifi dataset. List of hifi, e.g.: hf1.fq,hf2.fq')
# parser.add_argument('-c','--HIC', type=str, dest='Hi_C_dataset',
                    # help='List of Hi-C dataset(s), e.g.: hc1.fq,hc2.fq')
parser.add_argument('-t','--threads', type=int, dest='threads', default=4,
                    help='Number of threads, e.g.: 64')
parser.add_argument('-m','--ram', type=int, dest='ram', default=32,
                    help='Number of ram, minimum ram suggested: 32G')
parser.add_argument('-e','--extra_binner', type=str, dest='extra_binner',
                    help='Extra binner for binning: m: metabinner, v: vamb, l: lorbin; for instance: -e m, means BASALT will use metabinner for binning besides metabat2, maxbin2, and concoct')
parser.add_argument('-o','--out', type=str, dest='output_folder_name', default='Final_binset',
                    help='Name of the output folder. For binning, E.g. -o Anammox. BASALT would put those bins into folder Anammox_final_binset; for data feeding, e.g. -o Anammox; output files will under the folder of Anammox_data_feeded')
parser.add_argument('-q','--quality-check', type=str, dest='quality_check', default='checkm2', 
                    help='Chance checkm version, default: checkm2; you may use: \'-q checkm\' to specify checkm for quality check when running BASALT')
parser.add_argument('--min-cpn', type=int, dest='Min_completeness', default=35,
                    help='Min completeness of kept bins (default: 35)')
parser.add_argument('--max-ctn', type=int, dest='Max_contamination', default=20,
                    help='Max contamination of kept bins (default: 20)')
parser.add_argument('--mode', type=str, dest='running_mode', default='continue',
                    help='Start a new project (new) or continue to run (continue). e.g. --mode continue / --mode new')
parser.add_argument('--module', type=str, dest='functional_module', default='all',
                    help='Modules for binning. Four modules: 1. autobinning; 2. refinement; 3. reassembly; 4. all. Default will run all modules. But you could set the only perform modle. e.g. --module reassembly. In the module, ')
# parser.add_argument('--autopara', type=str, dest='autobining_parameters', default='more-sensitive',
#                     help='Three parameters to chose: 1. more-sensitive; 2. sensitive; 3. quick. Default: more-sensitive. e.g. --autopara sensitive')
parser.add_argument('--refinepara', type=str, dest='refinement_paramter', default='quick',
                    help='Two refinement parameters to chose: 1. deep; 2. quick. Default: quick. e.g. --refpara deep')
parser.add_argument('--sensitive', type=str, dest='binning_sensitive', default='sensitive',
                    help='Three parameters to chose: 1. quick; 2. sensitive; 3. more-sensitive. Default: sensitive. If you want to change the sensitive, use: e.g. --sensitive sensitive')
# parser.add_argument('--hybrid', type=str, dest='hybrid_reassembly', default='no',
#                     help='Use hybrid re-assembly. e.g. --hybrid y / --hybrid n; defalt no')
# parser.add_argument('--compression', type=str, dest='compression', default=0,
#                     help='Two refinement parameters to chose: 1. deep; 2. quick. Default: deep. e.g. --refpara quick')
# parser.add_argument('--hifi-only', action='store_true', default=False, dest='hifi',
#                     help='Only Hifi data')
# parser.add_argument('--ont-only', action='store_true', default=False, dest='ont',
#                     help='Only ont data')
# parser.add_argument('--hybrid', action='store_true', default=False, dest='ont',
#                     help='Hybrid data, including both HTS and ont/pb datasets')
parser.add_argument('-d','--data-feeding-folder', type=str, dest='data_feeding_folder',
                    help='List of folder name of extra binset(s) for data feeding, e.g.: -d binset1_folder_name,binset2_folder_name')
parser.add_argument('--binset-index', type=str, dest='extra_binset_start_index', default=500,
                    help='Optional parameter for data feeding. The start index of the extra binset, e.g.: -bi 5. BASALT already set a default index, but if you already had 4 assemblies, the binset start index could be 5')
# parser.add_argument('--only-refinement', action='store_true', dest='only_refinement',
#                     help='Only carry out refinement, e.g.: --only-refinement')
parser.add_argument('-r', '--refinement-binset', type=str, dest='refinement_binset', default='',
                    help='Specify binset folder name for refinement e.g.: -r Human_gut_microbime_MAGs')
parser.add_argument('-c', '--coverage-list', type=str, dest='coverage_list',
                    help='List of depth file for refinement. Coverage file(s) could be generated from data feeding modole. e.g.: -c Coverage_matrix_for_binning_1_assembly.fa.txt,Coverage_matrix_for_binning_2_assembly.fa.txt')
parser.add_argument('-b', '--binsets-list', type=str, dest='binsets_list',
                    help='List of binsets for de-replication. Binset depth file(s) could be generated from data feeding modole. e.g.: -b 1_assembly_BestBinsSet,2_assembly_BestBinsSet')

# args parsed inside main()
# args attributes consumed inside main()

# input normalisation moved into main()


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Banner
# ---------------------------------------------------------------------------

_BANNER = r"""
 ____    _    ____    _    _   _____
| __ )  / \  / ___|  / \  | | |_   _|
|  _ \ / _ \ \___ \ / _ \ | |   | |
| |_) / ___ \ ___) / ___ \| |___| |
|____/_/   \_\____/_/   \_\_____|_|
"""


def _print_banner():
    """Print the BASALT logo + version. CLI-only — never runs at ``import`` time."""
    from basalt import __version__
    line = _BANNER.rstrip('\n')
    tag = 'v' + __version__ + '  Metagenomic binning & refinement pipeline'
    print(line)
    print('  ' + tag)
    print()


def _fail(msg):
    sys.stderr.write('[BASALT] ERROR: ' + msg + '\n')
    sys.exit(1)


def _check_files_exist(paths, label):
    missing = [p for p in paths if p and not os.path.isfile(p)]
    if missing:
        _fail('{} file(s) not found: {}'.format(label, ', '.join(missing)))


def validate_inputs():
    # At least one input kind must be provided
    if not (assembly_list or binsets_list or refinement_binset or data_feeding_folder):
        _fail('no input provided: pass at least one of -a / -b / -r / -d')

    _check_files_exist(assembly_list, 'assembly')
    for pair in datasets.values():
        _check_files_exist(pair, 'short-read')
    _check_files_exist(lr_list, 'long-read')
    _check_files_exist(hifi_list, 'hifi-read')
    _check_files_exist(hic_list, 'Hi-C')
    _check_files_exist(coverage_list, 'coverage')

    for folder in (data_feeding_folder or []):
        if not os.path.isdir(folder):
            _fail('data-feeding folder not found: {}'.format(folder))
    for folder in (binsets_list or []):
        if not os.path.isdir(folder):
            _fail('binset folder not found: {}'.format(folder))
    if refinement_binset and not os.path.isdir(refinement_binset):
        _fail('refinement binset folder not found: {}'.format(refinement_binset))

    if not (1 <= num_threads <= 1024):
        _fail('--threads must be in [1, 1024]; got {}'.format(num_threads))
    if ram <= 0:
        _fail('--ram must be positive; got {}'.format(ram))
    if not (0 <= min_cpn <= 100):
        _fail('--min-cpn must be in [0, 100]; got {}'.format(min_cpn))
    if not (0 <= max_ctn <= 100):
        _fail('--max-ctn must be in [0, 100]; got {}'.format(max_ctn))

    if QC_software not in ('checkm', 'checkm2'):
        _fail('--quality-check must be "checkm" or "checkm2"; got {!r}'.format(QC_software))
    if continue_mode not in ('new', 'continue', 'last'):
        _fail('--mode must be "new" or "continue"; got {!r}'.format(continue_mode))
    if sensitivity not in ('quick', 'sensitive', 'more-sensitive'):
        _fail('--sensitive must be quick/sensitive/more-sensitive; got {!r}'.format(sensitivity))
    if refinement_paramter not in ('quick', 'deep'):
        _fail('--refinepara must be "quick" or "deep"; got {!r}'.format(refinement_paramter))
    if functional_module not in ('all', 'autobinning', 'refinement', 'reassembly'):
        _fail('--module must be all/autobinning/refinement/reassembly; got {!r}'.format(functional_module))

    # Output folder: verify parent directory is writable (BASALT suffixes
    # the stem itself, e.g. {stem}_final_binset, so don't pre-create it).
    out_parent = os.path.dirname(os.path.abspath(output_folder)) or '.'
    if not os.path.isdir(out_parent):
        _fail('output folder parent {!r} does not exist'.format(out_parent))
    if not os.access(out_parent, os.W_OK):
        _fail('output folder parent {!r} is not writable'.format(out_parent))


# validate_inputs() invoked inside main()


# config summary + continue-mode handling moved into main()

def _dispatch():
    """
    Entry point for the BASALT command-line interface.

    This function parses command-line arguments (already configured in the
    module-level ``parser``), normalises them into Python data structures,
    and dispatches to the appropriate BASALT workflow depending on:

    - quality control backend (CheckM2 vs. CheckM),
    - selected functional module (autobinning / refinement / reassembly / all),
    - whether data feeding or dereplication on existing binsets is requested.

    Returns
    -------
    None
        Side effects include running the end-to-end BASALT pipeline and
        writing results to the specified output folder.
    """
    global pwd, output_folder
    # Main execution: dispatch to the appropriate phase modules.
    # All phase modules accept a ``QC_software`` parameter ('checkm2' or
    # 'checkm') and configure the underlying step modules accordingly,
    # so we no longer need a top-level CheckM2 / CheckM branch here.

    if len(binsets_list) != 0:
        # Standalone S4 dereplication on externally-provided binsets.
        from basalt.steps.s4_multiple_assembly_comparator import (
            multiple_assembly_comparator_main,
            set_qc_backend as _set_s4_backend,
        )
        _set_s4_backend(QC_software)
        step = 'initial_drep'
        multiple_assembly_comparator_main(
            assembly_list, binsets_list, coverage_list, datasets, step, num_threads
        )

    elif refinement_binset != '':
        # Standalone S5 outlier removal on a refinement target.
        from basalt.steps.s5_outlier_remover_dl import (
            outlier_remover_main,
            set_qc_backend as _set_s5_backend,
        )
        _set_s5_backend(QC_software)
        outlier_remover_main(
            refinement_binset, coverage_list, datasets, assembly_list,
            pwd, num_threads, lr=lr_list, hifi_list=hifi_list,
        )

    elif len(assembly_list) != 0:
        # Full pipeline. Each phase module configures its own step backends
        # via the QC_software argument.
        if functional_module == 'autobinning' or functional_module == 'all':
            from basalt.module.autobinning import BASALT_main_autobinning
            BASALT_main_autobinning(
                assembly_list, datasets, num_threads, lr_list, hifi_list,
                hic_list, eb_list, ram, continue_mode, functional_module,
                sensitivity, refinement_paramter, max_ctn, min_cpn, pwd,
                QC_software, output_folder,
            )

        if functional_module == 'refinement' or functional_module == 'all':
            from basalt.module.refinement import BASALT_main_refinement
            BASALT_main_refinement(
                assembly_list, datasets, num_threads, lr_list, hifi_list,
                hic_list, eb_list, ram, continue_mode, functional_module,
                sensitivity, refinement_paramter, max_ctn, min_cpn, pwd,
                QC_software, output_folder,
            )

        if functional_module == 'reassembly' or functional_module == 'all':
            from basalt.module.reassembly import BASALT_main_re_assembly
            BASALT_main_re_assembly(
                assembly_list, datasets, num_threads, lr_list, hifi_list,
                hic_list, eb_list, ram, continue_mode, functional_module,
                sensitivity, refinement_paramter, max_ctn, min_cpn, pwd,
                QC_software, output_folder,
            )

        if len(data_feeding_folder) != 0:
            from basalt.module.datafeeding import data_feeding_main
            data_feeding_main(
                assembly_list, datasets, num_threads, data_feeding_folder,
                pwd, QC_software, output_folder, binsetindex, continue_mode,
            )

        from basalt.core.cleanup import cleanup
        cleanup(assembly_list)
        print('All accomplish!')

    else:
        # No assemblies — only data feeding requested. Call the low-level
        # Data_feeding routine directly (no S4/S5 post-processing without
        # an existing BASALT run to integrate into).
        if len(data_feeding_folder) != 0:
            from basalt.core.data_feeding import data_feeding
            if output_folder != 'Final_binset':
                output_folder = output_folder + '_data_feeded'
            else:
                output_folder = 'Data_feeded'
            data_feeding(
                data_feeding_folder, datasets, binsetindex,
                num_threads, output_folder, QC_software, pe='y',
            )


def main():
    """Console-script entry point. Parses ``sys.argv`` and dispatches."""
    _print_banner()
    global args
    global assemblies, sr_datasets, long_reads_list, hifi_list, QC_software
    global num_threads, ram, extrabinner, min_cpn, max_ctn, continue_mode
    global functional_module, refinement_paramter, output_folder, sensitivity
    global data_feeding_folder, binsetindex, refinement_binset, coverage_list
    global binsets_list, datasets, assembly_list, lr_list, hic_list, eb_list, pwd

    args = parser.parse_args()
    # Re-run all the setup performed at module level previously, now using
    # the freshly-parsed args:
    assemblies = args.assemblies
    sr_datasets = args.sr_datasets
    long_reads_list = args.long_datasets
    hifi_list = args.hifi_datasets
    QC_software = args.quality_check
    num_threads = args.threads
    ram = args.ram
    extrabinner = args.extra_binner
    min_cpn = args.Min_completeness
    max_ctn = args.Max_contamination
    continue_mode = args.running_mode
    functional_module = args.functional_module
    refinement_paramter = args.refinement_paramter
    output_folder = args.output_folder_name
    sensitivity = args.binning_sensitive
    data_feeding_folder = args.data_feeding_folder
    binsetindex = args.extra_binset_start_index
    refinement_binset = args.refinement_binset
    coverage_list = args.coverage_list
    binsets_list = args.binsets_list

    # Parse and normalise input lists.
    try:
        datasets_list = sr_datasets.split('/')
        datasets = {}
        for n, item in enumerate(datasets_list, start=1):
            pr = str(item).strip().split(',')
            datasets[str(n)] = [pr[0].strip(), pr[1].strip()]
    except Exception:
        datasets = {}
    try:
        assembly_list = assemblies.split(',')
    except Exception:
        assembly_list = []
    try:
        lr_list = [x.strip() for x in long_reads_list.split(',')]
    except Exception:
        lr_list = []
    try:
        hifi_list = [x.strip() for x in hifi_list.split(',')]
    except Exception:
        hifi_list = []
    hic_list = []  # Hi-C input not currently exposed
    try:
        eb_list = [x.strip() for x in extrabinner.split(',')]
    except Exception:
        eb_list = []
    try:
        data_feeding_folder = [x.strip() for x in data_feeding_folder.split(',')]
    except Exception:
        data_feeding_folder = []
    try:
        coverage_list = [x.strip() for x in coverage_list.split(',')]
    except Exception:
        coverage_list = []
    try:
        binsets_list = [x.strip() for x in binsets_list.split(',')]
    except Exception:
        binsets_list = []

    validate_inputs()

    pwd = os.getcwd()
    print('Processing assemblies:', str(assembly_list))
    print('Processing short-reads:', str(datasets))
    print('Processing long-reads:', str(lr_list))
    print('Processing hifi-reads:', str(hifi_list))
    print('Processing Hi-C reads:', str(hic_list))
    print('Output folder name will be:', str(output_folder))
    print('Process with extra binner:', str(eb_list))
    print('Quality check software:', str(QC_software))
    print('Binning sensitivity:', str(sensitivity))
    print('Processing with:', str(num_threads), 'threads')
    print('Processing with:', str(ram), 'G')
    print('Running status:', str(continue_mode))
    print('Binning module:', str(functional_module))
    print('Min completeness:', str(min_cpn))
    print('Max contamination:', str(max_ctn))
    print('Refinement parameter:', str(refinement_paramter))
    print('Extra binset(s) for data feeding:', str(data_feeding_folder))
    print('Refinement binset:', str(refinement_binset))
    print('List of coverage file(s):', str(coverage_list))
    print('Binset(s) list:', str(binsets_list))

    if continue_mode == 'continue':
        continue_mode = 'last'

    _dispatch()


if __name__ == '__main__':
    main()
