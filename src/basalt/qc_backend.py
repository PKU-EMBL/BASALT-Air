#!/usr/bin/env python

"""
Quality-check backend abstraction for BASALT.

BASALT historically shipped parallel copies of every S* module, one wired to
CheckM2 and one to CheckM, because the two tools differ in three places:

1. How they are invoked on the command line.
2. Where they write their per-bin result TSV.
3. The format of that TSV (CheckM2 is a flat table; CheckM is a dict-like
   text blob that embeds Python-dict-style metadata).

Using this module, a single S* implementation can delegate all three concerns
to a backend object and stop branching on ``QC_software`` at every call site.

Schema produced by ``parse_results``:

    {
        bin_id: {
            'Completeness':     float,   # 0-100
            'Contamination':    float,   # 0-100
            'Genome size':      int,     # bp
            'contig_size':      float,   # N50 (CheckM2) or Mean scaffold length (CheckM)
            'contig_size_key':  str,     # name of the metric: 'N50' or 'Mean scaffold length'
        },
        ...
    }

Call sites that historically read ``bin['N50']`` or ``bin['Mean scaffold length']``
should migrate to ``bin['contig_size']`` so they work under either backend.
"""

from __future__ import annotations

import os
from glob import glob


def _normalise_col(name):
    """Normalise a TSV column header for tolerant name matching.

    Lower-cases and strips spaces/underscores so that ``Genome_Size``,
    ``Genome size`` and ``genome_size`` all collapse to ``genomesize``.
    """
    return name.strip().lower().replace(' ', '').replace('_', '')


def get_backend(qc_software):
    """
    Factory: return the QC backend instance for the given software name.

    Parameters
    ----------
    qc_software : str
        Either ``'checkm2'`` or ``'checkm'``.
    """
    if qc_software == 'checkm2':
        return CheckM2Backend()
    if qc_software == 'checkm':
        return CheckMBackend()
    raise ValueError('unknown QC backend: {!r}'.format(qc_software))


class QCBackend:
    """Abstract base: concrete subclasses must implement all methods."""

    #: Identifier used in config/logs.
    name = 'abstract'
    #: Key name under which the per-contig/per-scaffold quality metric is
    #: commonly referred to in BASALT code. Schema always exposes this as
    #: ``contig_size``; this attribute names the underlying metric for logs.
    contig_size_key = 'contig_size'

    def build_cmd(self, threads, bin_folder, output_dir, ext='fa'):
        raise NotImplementedError

    def results_path(self, output_dir):
        """Path to the results TSV that ``parse_results`` will read."""
        raise NotImplementedError

    def parse_results(self, containing_folder):
        """
        Walk ``containing_folder`` and parse any results TSV found.

        Returns
        -------
        dict of bin_id -> metrics dict (see module docstring for schema).
        """
        raise NotImplementedError


class CheckM2Backend(QCBackend):
    """Backend for CheckM2 (``checkm2 predict`` → ``quality_report.tsv``)."""

    name = 'checkm2'
    contig_size_key = 'N50'
    results_filename = 'quality_report.tsv'

    def build_cmd(self, threads, bin_folder, output_dir, ext='fa'):
        return (
            'checkm2 predict -t {threads} -i {bin_folder} -x {ext} -o {output_dir}'
            .format(threads=threads, bin_folder=bin_folder, ext=ext, output_dir=output_dir)
        )

    def results_path(self, output_dir):
        return os.path.join(output_dir, self.results_filename)

    def parse_results(self, containing_folder):
        # Two layouts share this filename and both must be read:
        #   * native CheckM2 (14 cols): Name, Completeness, Contamination, ...,
        #     Contig_N50, Average_Gene_Length, Genome_Size, ...
        #   * BASALT's internal report (5 cols): Bin_ID, Genome_size,
        #     Completeness, Contamination, N50
        # We resolve columns by (normalised) header name and fall back to the
        # native CheckM2 positions when a header is absent/unrecognised, so an
        # older positional file still parses.
        bins = {}
        for root, _dirs, files in os.walk(containing_folder):
            for fname in files:
                if fname != self.results_filename:
                    continue
                with open(os.path.join(root, fname)) as fh:
                    lines = fh.readlines()
                if not lines:
                    continue

                header = [_normalise_col(h) for h in lines[0].rstrip('\n').split('\t')]
                col = {name: i for i, name in enumerate(header)}

                def _idx(candidates, fallback):
                    for c in candidates:
                        if c in col:
                            return col[c]
                    return fallback

                i_comp = _idx(('completeness',), 1)
                i_cont = _idx(('contamination',), 2)
                i_n50 = _idx(('contign50', 'n50'), 6)
                i_size = _idx(('genomesize',), 8)
                need = max(i_comp, i_cont, i_n50, i_size)

                for line in lines[1:]:
                    fields = line.rstrip('\n').split('\t')
                    if len(fields) <= need:
                        continue
                    bin_id = fields[0]
                    if '_genomes.0' in bin_id:
                        continue
                    try:
                        completeness = float(fields[i_comp])
                        contamination = float(fields[i_cont])
                        n50 = float(fields[i_n50])
                        genome_size = int(float(fields[i_size]))
                    except (ValueError, IndexError):
                        continue
                    bins[bin_id] = {
                        'Completeness': completeness,
                        'Contamination': contamination,
                        'Genome size': genome_size,
                        'contig_size': n50,
                        'contig_size_key': 'N50',
                        # Legacy aliases for call sites still reading the raw key:
                        'N50': n50,
                    }
        return bins


class CheckMBackend(QCBackend):
    """Backend for CheckM (``checkm lineage_wf`` → ``storage/bin_stats_ext.tsv``)."""

    name = 'checkm'
    contig_size_key = 'Mean scaffold length'
    results_filename = 'bin_stats_ext.tsv'

    def build_cmd(self, threads, bin_folder, output_dir, ext='fa'):
        return (
            'checkm lineage_wf -t {threads} -x {ext} {bin_folder} {output_dir}'
            .format(threads=threads, ext=ext, bin_folder=bin_folder, output_dir=output_dir)
        )

    def results_path(self, output_dir):
        return os.path.join(output_dir, 'storage', self.results_filename)

    def parse_results(self, containing_folder):
        """
        CheckM writes a quirky TSV where column 2 is a Python-dict-ish string
        like ``{'marker lineage': '...', 'Completeness': 95.1, ...}``. Parse
        it by ad-hoc string splitting (matches legacy behaviour) rather than
        eval for safety.
        """
        bins = {}
        for root, _dirs, files in os.walk(containing_folder):
            for fname in files:
                if fname != self.results_filename:
                    continue
                with open(os.path.join(root, fname)) as fh:
                    for line in fh:
                        bin_id = line.strip().split('\t')[0]
                        if '_genomes.0' in bin_id:
                            continue
                        metrics = {}
                        try:
                            metrics['marker lineage'] = line.split("'marker lineage': '")[1].split("'")[0]
                        except IndexError:
                            metrics['marker lineage'] = 'root'
                        try:
                            metrics['Completeness'] = float(
                                line.split("'Completeness': ")[1].split(', ')[0]
                            )
                        except (IndexError, ValueError):
                            metrics['Completeness'] = 0.0
                        try:
                            metrics['Genome size'] = int(float(
                                line.split("'Genome size':")[1].split(', ')[0]
                                    .replace("'", '').replace('"', '').replace('}', '').strip()
                            ))
                        except (IndexError, ValueError):
                            metrics['Genome size'] = 0
                        try:
                            metrics['Contamination'] = float(
                                line.split("Contamination': ")[1].split('}')[0].split(',')[0]
                            )
                        except (IndexError, ValueError):
                            metrics['Contamination'] = 0.0
                        try:
                            msl_raw = line.split("Mean scaffold length':")[1]
                            # Works for both ``...', key: ...`` and ``...'}}`` terminators.
                            msl_str = msl_raw.split(',')[0].split('}')[0].strip()
                            msl = float(msl_str)
                        except (IndexError, ValueError):
                            msl = 0.0
                        metrics['Mean scaffold length'] = msl
                        metrics['contig_size'] = msl
                        metrics['contig_size_key'] = 'Mean scaffold length'
                        bins[bin_id] = metrics
        return bins


if __name__ == '__main__':
    # Minimal sanity check: factory returns the expected classes.
    assert isinstance(get_backend('checkm2'), CheckM2Backend)
    assert isinstance(get_backend('checkm'), CheckMBackend)
    print('qc_backend.py sanity check OK')
