# BASALT-Air v1.0.0

**Binning Across a Series of Assemblies Toolkit, Air version.**

BASALT-Air is a lightweight, modular metagenomic binning and refinement pipeline for generating high-quality metagenome-assembled genomes (MAGs) from short reads, long reads, and hybrid assemblies.

[![Nature Communications](https://img.shields.io/badge/Nature%20Communications-2024-blue)](https://doi.org/10.1038/s41467-024-46539-7)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Highlights

- Absolute input path support, so runs do not need to start inside data directories.
- Multi-assembly binning and dereplication.
- Deep-learning-based bin refinement.
- Short-read, ONT/PacBio long-read, and PacBio HiFi support.
- Modular execution: autobinning, refinement, reassembly, and datafeeding.
- Reproducible run manifests and timestamped logs.

The command-line entry points installed by this package are lowercase:

```bash
basalt
basalt_models_download
```

The internal Python package remains `basalt` for compatibility with the existing codebase.

## Workflow

<p align="center">
  <img src="fig/workflow.png" width="75%" alt="BASALT-Air workflow">
</p>

BASALT-Air consists of four modules:

| Module | Steps | Description |
| --- | --- | --- |
| `autobinning` | S1-S4 | Read mapping, binning, within-group dereplication, and multi-assembly dereplication |
| `refinement` | S5-S7 | Deep-learning outlier removal, contig retrieval, and optional polishing |
| `reassembly` | S8-S10 | OLC elongation, short-read/hybrid reassembly, and post-reassembly dereplication |
| `datafeeding` | - | Integrate external binsets and rerun downstream comparison/refinement |

Default `--module all` runs autobinning, refinement, and reassembly sequentially.

## Installation

BASALT-Air requires Linux and Python 3.12. The recommended installation path is [pixi](https://pixi.sh), which installs Python dependencies and the external bioinformatics tools listed in `pixi.toml`.

### 1. Install Pixi

```bash
curl -fsSL https://pixi.sh/install.sh | sh
```

### 2. Clone The Repository

```bash
git clone https://github.com/PKU-EMBL/BASALT-Air.git
cd BASALT-Air
```

### 3. Configure Local Paths

Edit `pixi.toml` before the first run:

```toml
[activation.env]
BASALT_WEIGHT = "/your/path/to/basalt_weights"
CHECKM2DB     = "/your/path/to/checkm2db/CheckM2_database/uniref100.KO.1.dmnd"
```

If your system needs a different CUDA compatibility target, update:

```toml
[system-requirements]
cuda = "12"
```

### 4. Install The Environment

```bash
pixi install
```

### 5. Download Model Weights And Databases

BASALT-Air model weights are available from:

- Hugging Face: <https://huggingface.co/PKU-EMBL/BASALT_WEIGHT>
- Google Drive: <https://drive.google.com/drive/folders/1d0e_2FpYRBAZLwKXl8fA-yDK4b5PBA_E>
- Baidu Netdisk: <https://pan.baidu.com/s/1ouKqabxHYr1GmvpquQCzqw?pwd=embl> (提取码: `embl`)

CheckM2 database and demo data are available from the same Google Drive and Baidu Netdisk links.

Expected local paths:

- `basalt_weights/` should match `BASALT_WEIGHT`.
- `checkm2db/CheckM2_database/uniref100.KO.1.dmnd` should match `CHECKM2DB`.
- `checkmdb/` is optional and only needed for legacy CheckM runs.

You can download model weights with either Hugging Face CLI:

```bash
pip install huggingface_hub
huggingface-cli download PKU-EMBL/BASALT_WEIGHT --local-dir /your/path/to/basalt_weights
```

or the pixi task:

```bash
pixi run download-weights
```

Download the CheckM2 database with:

```bash
pixi run checkm2-db
```

## Verify Installation

```bash
pixi shell
basalt --version
basalt --check-deps
```

You can also run the predefined tasks without entering the shell:

```bash
pixi run version
pixi run check-deps
pixi run sanity
```

## Quick Start

Single assembly with paired-end reads:

```bash
basalt -a assembly.fa -s r1.fq,r2.fq -t 32 -m 128
```

Multiple assemblies and datasets:

```bash
basalt -a as1.fa,as2.fa,as3.fa \
       -s d1_r1.fq,d1_r2.fq/d2_r1.fq,d2_r2.fq/d3_r1.fq,d3_r2.fq \
       -t 60 -m 250
```

Hybrid assembly with short reads, ONT/PacBio long reads, and HiFi reads:

```bash
basalt -a assembly.fa \
       -s sr_r1.fq,sr_r2.fq \
       -l ont.fq \
       -hf hifi.fq \
       -t 60 -m 250
```

Run from any directory with absolute paths:

```bash
basalt \
    -a /path/to/data/assembly.fa \
    -s /path/to/data/sample1.R1.fq,/path/to/data/sample1.R2.fq \
    -l /path/to/data/sample1.nanopore.fq \
    -t 64 -m 128 \
    -o my_project \
    --workdir /scratch/work \
    --outdir /results/output
```

For multiple absolute-path paired-end datasets, use `;` between read pairs:

```bash
basalt \
    -a /data/as1.fa,/data/as2.fa \
    -s /data/s1_R1.fq,/data/s1_R2.fq;/data/s2_R1.fq,/data/s2_R2.fq \
    -t 64 -m 128
```

## Common Options

| Flag | Description | Default |
| --- | --- | --- |
| `-a`, `--assemblies` | Assembly FASTA files, comma-separated | - |
| `-s`, `--shortreads` | Paired-end short reads | - |
| `-l`, `--longreads` | ONT/PacBio CLR long reads | - |
| `-hf`, `--hifi` | PacBio HiFi reads | - |
| `-t`, `--threads` | Number of threads | `4` |
| `-m`, `--ram` | RAM in GB | `32` |
| `-q`, `--quality-check` | `checkm2` or `checkm` | `checkm2` |
| `--module` | `all`, `autobinning`, `refinement`, or `reassembly` | `all` |
| `--min-cpn` | Minimum completeness percentage | `35` |
| `--max-ctn` | Maximum contamination percentage | `20` |
| `-o`, `--out` | Output folder name | `Final_binset` |
| `--workdir` | Directory for intermediate files | current directory |
| `--outdir` | Directory for final output | same as `--workdir` |

Run `basalt --help` for the full option list.

## Demo Dataset

Demo data are available from:

- Google Drive: <https://drive.google.com/drive/folders/1d0e_2FpYRBAZLwKXl8fA-yDK4b5PBA_E>
- Baidu Netdisk: <https://pan.baidu.com/s/1ouKqabxHYr1GmvpquQCzqw?pwd=embl> (提取码: `embl`)

Example:

```bash
tar -xzf BASALT_demo.tar.gz

basalt \
    -a /path/to/BASALT_demo/Data/assembly.fa \
    -s /path/to/BASALT_demo/Data/sample1.R1.fq,/path/to/BASALT_demo/Data/sample1.R2.fq \
    -l /path/to/BASALT_demo/Data/sample1.nanopore.fq \
    -t 64 -m 128 \
    -o demo_run \
    --workdir /scratch/demo_work \
    --outdir /results/demo_output
```

## Output

Results are written to `<output_folder>/`, defaulting to `Final_binset/`.

Key output files include:

- `*.fa`: final dereplicated MAGs.
- `OLC_quality_report.tsv`: completeness, contamination, N50, and related metrics.
- `BASALT_run_manifest.json`: reproducibility metadata.
- `Basalt_log.txt`: timestamped pipeline log.

## Citation

If you use BASALT-Air in your research, cite the original BASALT publication:

```bibtex
@article{qiu2024basalt,
  title   = {BASALT refines binning from metagenomic data and increases
             resolution of genome-resolved metagenomic analysis},
  author  = {Qiu, Zhiguang and Yuan, Li and Lian, Chun-Ang and Lin, Bin
             and Chen, Jie and Mu, Rong and Qiao, Xuejiao and Zhang, Liyu
             and Xu, Zheng and Fan, Lu and others},
  journal = {Nature Communications},
  volume  = {15},
  number  = {1},
  pages   = {2179},
  year    = {2024},
  doi     = {10.1038/s41467-024-46539-7}
}
```

## License

MIT License. See [LICENSE](LICENSE).

## Contact

For issues and questions, open an issue at <https://github.com/PKU-EMBL/BASALT-Air/issues> or contact <yuke.sz@pku.edu.cn> and <zrjiang25@stu.pku.edu.cn>.
