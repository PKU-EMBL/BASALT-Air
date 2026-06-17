#!/usr/bin/python
# -*- coding: UTF-8 -*-

"""
Command-line interface for the BASALT-Air metagenomic binning pipeline.

This script parses user arguments and dispatches to the corresponding
pipeline modules (autobinning, refinement, reassembly, data feeding,
and dereplication) depending on the selected mode and options.
"""

import time
import sys
import os
import argparse
import shutil
import warnings
from glob import glob

from basalt.logger import setup_logger, get_logger, format_elapsed

try:
    from sklearn.exceptions import EfficiencyWarning
    warnings.filterwarnings("ignore", category=EfficiencyWarning)
except ImportError:
    pass


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
_USAGE_EXAMPLES = """
examples
--------

  # Single-sample binning, paired-end short reads, CheckM2 QC (default):
  basalt -a assembly.fa -s r1.fq,r2.fq -t 32 -m 64 -o my_run

  # Multi-sample with absolute paths (use ';' to separate read pairs):
  basalt -a /data/asm1.fa,/data/asm2.fa \\
         -s /data/r1_1.fq,/data/r1_2.fq;/data/r2_1.fq,/data/r2_2.fq \\
         --workdir /scratch/work --outdir /results -o multi_sample

  # Hybrid assembly + long reads:
  basalt -a asm.fa -s r1.fq,r2.fq -l ont.fq -hf hifi.fq -t 64 -m 128

  # Re-run only the refinement module on existing binsets:
  basalt -r my_binset -c coverage.txt --module refinement -t 32

  # Verify external tools are installed before launching:
  basalt --check-deps

  # Preview the resolved configuration without running anything:
  basalt -a asm.fa -s r1.fq,r2.fq --dry-run
"""

parser = argparse.ArgumentParser(
    prog='basalt',
    description='BASALT-Air: metagenomic binning, refinement, and reassembly pipeline.',
    epilog=_USAGE_EXAMPLES,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
parser.add_argument('-V', '--version', action='version',
                    version='%(prog)s (BASALT-Air) ' + __import__('basalt').__version__)
parser.add_argument('--check-deps', action='store_true', dest='check_deps',
                    help='Check that required external tools (bowtie2, samtools, minimap2, metabat2, checkm2/checkm, blastn, ...) are on PATH and exit. Useful before submitting a long job.')
parser.add_argument('--dry-run', action='store_true', dest='dry_run',
                    help='Resolve all inputs/parameters, write the run manifest, then exit before running the pipeline. Useful to validate a command-line.')
parser.add_argument('-a', '--assemblies', type=str, dest='assemblies',
                    help='List of assemblies, e.g.: as1.fa,as2.fa')
parser.add_argument('-s', '--shortreads', type=str, dest='sr_datasets',
                    help='List of paired-end reads. Each pair is "mate1,mate2". Multiple pairs are separated by "/" or ";" (use ";" if any path is absolute, since "/" collides with directory separators). e.g.: r1_1.fq,r1_2.fq/r2_1.fq,r2_2.fq  or  /abs/r1_1.fq,/abs/r1_2.fq;/abs/r2_1.fq,/abs/r2_2.fq')
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
                    help='Extra binner(s) to run alongside metabat2, maxbin2 and concoct. v: vamb, l: lorbin. e.g. -e v or -e v,l')
parser.add_argument('-o','--out', type=str, dest='output_folder_name', default='Final_binset',
                    help='Name of the output folder. For binning, E.g. -o Anammox. BASALT-Air would put those bins into folder Anammox_final_binset; for data feeding, e.g. -o Anammox; output files will under the folder of Anammox_data_feeded')
parser.add_argument('-w','--workdir', type=str, dest='workdir', default=None,
                    help='Directory for intermediate process files. Default: current directory. BASALT-Air will create it if missing and chdir into it before running.')
parser.add_argument('--outdir', type=str, dest='outdir', default=None,
                    help='Directory where the final output folder is moved to once the pipeline finishes. Default: same as --workdir (output stays alongside intermediates).')
parser.add_argument('-q','--quality-check', type=str, dest='quality_check', default='checkm2', 
                    help='Chance checkm version, default: checkm2; you may use: \'-q checkm\' to specify checkm for quality check when running BASALT-Air')
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
                    help='Optional parameter for data feeding. The start index of the extra binset, e.g.: -bi 5. BASALT-Air already set a default index, but if you already had 4 assemblies, the binset start index could be 5')
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
 ____    _    ____    _    _   _____        _    ___ ____
| __ )  / \  / ___|  / \  | | |_   _|      / \  |_ _|  _ \
|  _ \ / _ \ \___ \ / _ \ | |   | | _____ / _ \  | || |_) |
| |_) / ___ \ ___) / ___ \| |___| ||_____/ ___ \ | ||  _ <
|____/_/   \_\____/_/   \_\_____|_|     /_/   \_\___|_| \_\
"""


def _print_banner():
    """Print the BASALT-Air logo + version. CLI-only — never runs at ``import`` time."""
    from basalt import __version__
    line = _BANNER.rstrip('\n')
    tag = 'BASALT-Air v' + __version__ + '  Metagenomic binning & refinement pipeline'
    print(line)
    print('  ' + tag)
    print()


def _fail(msg):
    sys.stderr.write('[BASALT-Air] ERROR: ' + msg + '\n')
    sys.exit(1)


# ---------------------------------------------------------------------------
# Dependency check
# ---------------------------------------------------------------------------

# Tools that are essentially always invoked.
_DEPS_CORE = [
    'bowtie2', 'bowtie2-build', 'samtools', 'minimap2',
    'blastn', 'makeblastdb',
    'metabat2', 'run_MaxBin.pl', 'concoct',
    'jgi_summarize_bam_contig_depths',
    'spades.py', 'pilon', 'perl',
]
_DEPS_QC = {'checkm2': ['checkm2'], 'checkm': ['checkm']}
_DEPS_EXTRA_BINNER = {'v': ['vamb'], 'l': ['LorBin']}


def _which(name):
    return shutil.which(name)


def _check_dependencies(qc_software='checkm2', extra_binners=()):
    """Return (missing, present) lists for the active tool set."""
    needed = list(_DEPS_CORE)
    needed += _DEPS_QC.get(qc_software, [])
    for code in extra_binners or ():
        needed += _DEPS_EXTRA_BINNER.get(code, [])
    seen, ordered = set(), []
    for t in needed:
        if t not in seen:
            seen.add(t)
            ordered.append(t)
    missing = [t for t in ordered if _which(t) is None]
    present = [t for t in ordered if _which(t) is not None]
    return missing, present


def _print_dep_report(qc_software, extra_binners):
    missing, present = _check_dependencies(qc_software, extra_binners)
    print('Dependency check (QC backend: {}, extra binners: {}):'.format(
        qc_software, list(extra_binners) or 'none'))
    for t in present:
        print('  [ OK ] {} -> {}'.format(t, _which(t)))
    for t in missing:
        print('  [MISS] {}'.format(t))
    return missing


# ---------------------------------------------------------------------------
# Run manifest (reproducibility)
# ---------------------------------------------------------------------------

def _build_manifest(args_namespace, resolved):
    """Capture argv, version, host, params and start time as a dict."""
    import getpass
    import platform
    import socket
    from basalt import __version__ as _v

    return {
        'basalt_version': _v,
        'argv': list(sys.argv),
        'cwd_at_invocation': os.getcwd(),
        'workdir': resolved.get('workdir'),
        'outdir': resolved.get('outdir'),
        'parameters': vars(args_namespace),
        'normalised_inputs': {
            'assemblies': resolved.get('assemblies'),
            'short_reads': resolved.get('short_reads'),
            'long_reads': resolved.get('long_reads'),
            'hifi_reads': resolved.get('hifi_reads'),
            'extra_binner': resolved.get('extra_binner'),
            'coverage_list': resolved.get('coverage_list'),
            'binsets_list': resolved.get('binsets_list'),
            'data_feeding_folder': resolved.get('data_feeding_folder'),
            'refinement_binset': resolved.get('refinement_binset'),
        },
        'host': {
            'hostname': socket.gethostname(),
            'platform': platform.platform(),
            'python': sys.version.split()[0],
            'user': getpass.getuser(),
            'pid': os.getpid(),
        },
        'started_at': time.strftime('%Y-%m-%dT%H:%M:%S%z'),
        'started_at_epoch': time.time(),
    }


def _write_manifest(manifest, workdir):
    import json
    path = os.path.join(workdir, 'BASALT_run_manifest.json')
    try:
        with open(path, 'w') as fh:
            json.dump(manifest, fh, indent=2, default=str)
    except OSError as e:
        sys.stderr.write(
            '[BASALT-Air] WARNING: could not write run manifest {!r}: {}\n'.format(path, e))
    return path


def _finalise_manifest(path, status, elapsed_seconds, error=None):
    """Append end_status / duration to the manifest written at start."""
    import json
    try:
        with open(path) as fh:
            data = json.load(fh)
    except (OSError, ValueError):
        return
    data['ended_at'] = time.strftime('%Y-%m-%dT%H:%M:%S%z')
    data['duration_seconds'] = round(float(elapsed_seconds), 3)
    data['status'] = status
    if error is not None:
        data['error'] = str(error)
    try:
        with open(path, 'w') as fh:
            json.dump(data, fh, indent=2, default=str)
    except OSError:
        pass


def _check_files_exist(paths, label):
    missing = [p for p in paths if p and not os.path.isfile(p)]
    if missing:
        _fail('{} file(s) not found: {}'.format(label, ', '.join(missing)))


def _parse_paired_reads(value):
    """Parse the -s/--shortreads value into ``{idx: [mate1, mate2]}``.

    Supports three forms:
      * legacy '/'-separated pair groups, e.g. ``a,b/c,d``;
      * ';'-separated pair groups (use this when paths are absolute, since
        '/' collides with the filesystem separator), e.g. ``/x/a,/x/b;/x/c,/x/d``;
      * a flat comma list — auto-paired by twos — used when '/' would
        collide with absolute paths and no ';' was given.
    Empty or missing input yields an empty dict (callers may dispatch a
    flow that needs no short reads).
    """
    if not value:
        return {}
    text = str(value).strip()
    if not text:
        return {}

    if ';' in text:
        groups = [g.strip() for g in text.split(';') if g.strip()]
    elif '/' in text and not text.lstrip().startswith('/'):
        # Legacy form: "a,b/c,d". Only safe when no token is absolute.
        groups = [g.strip() for g in text.split('/') if g.strip()]
    else:
        # Either no separator at all (single pair "a,b") or the input
        # contains absolute paths whose '/' would corrupt a '/'-split.
        # Treat as a flat comma list and pair adjacent items.
        tokens = [t.strip() for t in text.split(',') if t.strip()]
        if len(tokens) % 2 != 0:
            _fail('--shortreads must list an even number of files when '
                  'pairs are not delimited by "/" or ";"; got {} file(s)'
                  .format(len(tokens)))
        groups = [','.join(tokens[i:i + 2]) for i in range(0, len(tokens), 2)]

    pairs = {}
    for n, item in enumerate(groups, start=1):
        files = [x.strip() for x in item.split(',') if x.strip()]
        if len(files) != 2:
            _fail('--shortreads pair {} must have exactly 2 files; got: {!r}'
                  .format(n, item))
        pairs[str(n)] = files
    return pairs


def _link_into_cwd(path):
    """Make ``path`` reachable as a bare filename in the current directory.

    Downstream BASALT modules reference inputs by basename (e.g.
    ``pwd + '/' + name``, ``unzip name``, ``os.chdir(name)``), so absolute
    or sub-directory paths are normalised here by symlinking them into cwd
    and returning the link's basename. Bare filenames pass through unchanged.
    """
    if not path:
        return path
    has_sep = (os.sep in path) or (os.altsep is not None and os.altsep in path)
    if not has_sep:
        return path
    abs_path = os.path.abspath(os.path.expanduser(path))
    name = os.path.basename(abs_path.rstrip(os.sep))
    if not name:
        return path
    cwd = os.getcwd()
    link = os.path.join(cwd, name)
    if os.path.realpath(link) == abs_path:
        return name
    if os.path.lexists(link):
        if os.path.islink(link):
            try:
                os.unlink(link)
            except OSError as e:
                _fail('cannot replace stale symlink {!r}: {}'.format(link, e))
        else:
            sys.stderr.write(
                '[BASALT-Air] WARNING: cwd already contains {!r}; '
                'using it instead of {!r}\n'.format(name, abs_path))
            return name
    try:
        os.symlink(abs_path, link)
    except OSError as e:
        _fail('failed to symlink {!r} -> {!r}: {}'.format(link, abs_path, e))
    return name


def _normalize_paths_to_cwd():
    """Symlink every path-bearing input into cwd so the rest of the pipeline,
    which assumes bare filenames in the working directory, keeps working
    when the user passes absolute or sub-directory paths."""
    global assembly_list, lr_list, hifi_list, hic_list
    global datasets, coverage_list, binsets_list
    global data_feeding_folder, refinement_binset

    assembly_list = [_link_into_cwd(p) for p in assembly_list]
    lr_list = [_link_into_cwd(p) for p in lr_list]
    hifi_list = [_link_into_cwd(p) for p in hifi_list]
    hic_list = [_link_into_cwd(p) for p in hic_list]
    coverage_list = [_link_into_cwd(p) for p in coverage_list]
    for k, pair in list(datasets.items()):
        datasets[k] = [_link_into_cwd(pair[0]), _link_into_cwd(pair[1])]

    data_feeding_folder = [_link_into_cwd(p) for p in data_feeding_folder]
    binsets_list = [_link_into_cwd(p) for p in binsets_list]
    if refinement_binset:
        refinement_binset = _link_into_cwd(refinement_binset)


def _move_output_to_outdir():
    """Move the produced output folder from workdir to outdir.

    The pipeline writes its final binset to ``<workdir>/<output_folder>``.
    When ``--outdir`` was supplied (and differs from workdir), relocate that
    folder so users can keep bulky intermediates separate from results.
    """
    src = os.path.join(workdir, output_folder)
    if not os.path.isdir(src):
        sys.stderr.write(
            '[BASALT-Air] WARNING: expected output folder {!r} not found in '
            'workdir; nothing to move to --outdir\n'.format(output_folder))
        return
    dst = os.path.join(outdir, output_folder)
    if os.path.exists(dst):
        sys.stderr.write(
            '[BASALT-Air] WARNING: --outdir already contains {!r}; left output '
            'in workdir at {!r}\n'.format(output_folder, src))
        return
    try:
        shutil.move(src, dst)
    except (OSError, shutil.Error) as e:
        sys.stderr.write(
            '[BASALT-Air] WARNING: could not move {!r} to {!r}: {}; left in '
            'workdir\n'.format(src, dst, e))
        return
    print('Final output moved to', dst)


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

    # -e/--extra_binner: only 'v' (VAMB) and 'l' (LorBin) are supported.
    # 'm' (MetaBinner) was removed in v1.2.x.
    for code in (eb_list or []):
        if code == 'm':
            _fail("MetaBinner ('-e m') is no longer supported. Use '-e v' (VAMB) "
                  "and/or '-e l' (LorBin) instead, or omit -e entirely.")
        if code not in ('v', 'l'):
            _fail("--extra_binner code must be 'v' (vamb) or 'l' (lorbin); got {!r}".format(code))


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
        Side effects include running the end-to-end BASALT-Air pipeline and
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
    global workdir, outdir

    args = parser.parse_args()

    # --check-deps: verify external tools are on PATH and exit (0 if all
    # found, 1 if any missing). Honours the QC backend and -e selection so
    # users only check what's actually needed for their planned run.
    if getattr(args, 'check_deps', False):
        eb_codes = []
        if args.extra_binner:
            eb_codes = [c.strip() for c in args.extra_binner.split(',') if c.strip()]
        missing = _print_dep_report(args.quality_check, eb_codes)
        sys.exit(1 if missing else 0)

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
    workdir = args.workdir
    outdir = args.outdir
    sensitivity = args.binning_sensitive
    data_feeding_folder = args.data_feeding_folder
    binsetindex = args.extra_binset_start_index
    refinement_binset = args.refinement_binset
    coverage_list = args.coverage_list
    binsets_list = args.binsets_list

    # Parse and normalise input lists.
    datasets = _parse_paired_reads(sr_datasets)
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

    # Allow path-bearing -o for backward compat: split into outdir + basename.
    if output_folder:
        _od, _bn = os.path.split(output_folder)
        if _od:
            if outdir is None:
                outdir = _od
            output_folder = _bn or 'Final_binset'

    # Resolve workdir (where intermediate process files live) and chdir into
    # it so the rest of the pipeline — which assumes cwd is the work dir —
    # keeps working unchanged.
    if workdir:
        workdir = os.path.abspath(os.path.expanduser(workdir))
        try:
            os.makedirs(workdir, exist_ok=True)
        except OSError as e:
            _fail('cannot create --workdir {!r}: {}'.format(workdir, e))
        os.chdir(workdir)
    else:
        workdir = os.getcwd()
    if not os.access(workdir, os.W_OK):
        _fail('working directory {!r} is not writable'.format(workdir))

    # Resolve outdir (where the final output folder is moved to). Defaults
    # to workdir, in which case no post-pipeline move happens.
    if outdir:
        outdir = os.path.abspath(os.path.expanduser(outdir))
        try:
            os.makedirs(outdir, exist_ok=True)
        except OSError as e:
            _fail('cannot create --outdir {!r}: {}'.format(outdir, e))
        if not os.access(outdir, os.W_OK):
            _fail('--outdir {!r} is not writable'.format(outdir))
    else:
        outdir = workdir

    pwd = workdir
    _normalize_paths_to_cwd()

    # Configure the unified logger now that workdir is final. This makes
    # `Basalt_log.txt` land in the workdir alongside the existing manual writes.
    setup_logger(os.path.join(workdir, 'Basalt_log.txt'))
    log = get_logger()

    log.info('Processing assemblies: %s', assembly_list)
    log.info('Processing short-reads: %s', datasets)
    log.info('Processing long-reads: %s', lr_list)
    log.info('Processing hifi-reads: %s', hifi_list)
    log.info('Processing Hi-C reads: %s', hic_list)
    log.info('Output folder name will be: %s', output_folder)
    log.info('Working directory (intermediate files): %s', workdir)
    log.info('Final output directory: %s', outdir)
    log.info('Process with extra binner: %s', eb_list)
    log.info('Quality check software: %s', QC_software)
    log.info('Binning sensitivity: %s', sensitivity)
    log.info('Processing with: %s threads', num_threads)
    log.info('Processing with: %s G', ram)
    log.info('Running status: %s', continue_mode)
    log.info('Binning module: %s', functional_module)
    log.info('Min completeness: %s', min_cpn)
    log.info('Max contamination: %s', max_ctn)
    log.info('Refinement parameter: %s', refinement_paramter)
    log.info('Extra binset(s) for data feeding: %s', data_feeding_folder)
    log.info('Refinement binset: %s', refinement_binset)
    log.info('List of coverage file(s): %s', coverage_list)
    log.info('Binset(s) list: %s', binsets_list)

    if continue_mode == 'continue':
        continue_mode = 'last'

    # Auto dependency check before any work — fail fast with a clear list of
    # missing executables instead of dying mid-pipeline with a cryptic error.
    missing_deps, _ = _check_dependencies(QC_software, eb_list)
    if missing_deps:
        log.warning(
            'Missing external tools on PATH: %s. Run "basalt --check-deps" '
            'for the full report. Continuing — failures will surface when '
            'BASALT-Air shells out to the missing tool.',
            ', '.join(missing_deps),
        )

    # Reproducibility: capture argv, version, host, parameters, timing.
    manifest = _build_manifest(args, {
        'workdir': workdir, 'outdir': outdir,
        'assemblies': assembly_list, 'short_reads': datasets,
        'long_reads': lr_list, 'hifi_reads': hifi_list,
        'extra_binner': eb_list, 'coverage_list': coverage_list,
        'binsets_list': binsets_list,
        'data_feeding_folder': data_feeding_folder,
        'refinement_binset': refinement_binset,
    })
    manifest_path = _write_manifest(manifest, workdir)
    log.info('Run manifest written: %s', manifest_path)

    if getattr(args, 'dry_run', False):
        log.info('--dry-run: configuration validated; exiting before _dispatch()')
        return

    t_start = time.monotonic()
    log.info('BASALT-Air pipeline started')
    try:
        _dispatch()

        if outdir != workdir:
            _move_output_to_outdir()
    except BaseException as exc:
        elapsed = time.monotonic() - t_start
        log.error(
            'BASALT-Air pipeline aborted after %s (%s: %s)',
            format_elapsed(elapsed), type(exc).__name__, exc,
        )
        _finalise_manifest(manifest_path, 'aborted', elapsed, error=exc)
        raise
    else:
        elapsed = time.monotonic() - t_start
        log.info('BASALT-Air pipeline finished in %s (%.1f s)', format_elapsed(elapsed), elapsed)
        _finalise_manifest(manifest_path, 'success', elapsed)


if __name__ == '__main__':
    main()
