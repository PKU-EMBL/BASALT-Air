# BASALT-Air — Binning Across a Series of Assemblies Toolkit Air Version

```
 ____    _    ____    _    _   _____
| __ )  / \  / ___|  / \  | | |_   _|
|  _ \ / _ \ \___ \ / _ \ | |   | |
| |_) / ___ \ ___) / ___ \| |___| |
|____/_/   \_\____/_/   \_\_____|_|
       Metagenomic binning & refinement pipeline
```

BASALT-Air is a lighted, versatile, modular pipeline that produces high-quality
metagenome-assembled genomes (MAGs) from short reads, long reads, hybrid
assemblies, and pre-existing binsets. It runs autobinning, multi-assembly
dereplication, deep-learning–based outlier removal, contig retrieval, OLC
elongation, and reassembly in a single coherent workflow.

📜 Published in *Nature Communications* (2024) —
[doi:10.1038/s41467-024-46539-7](https://doi.org/10.1038/s41467-024-46539-7).

---

## 📣 News

* **[2026/05/18]** 🔧 **Configuration improvements for easier setup.**
  * **User-configurable paths:** `pixi.toml` now uses placeholder paths (`/path/to/...`) instead of hardcoded personal paths. Users must edit `BASALT_WEIGHT` and `CHECKM2DB` in `pixi.toml` before running `pixi install`. See the updated installation section for detailed instructions.
  * **CUDA version flexibility:** Added clear comments in `pixi.toml` explaining how to adjust the CUDA version (line 13) to match your system (`"11"`, `"12"`, or `"13"`).
  * **libopenblas dependency:** Explicitly added `libopenblas` to dependencies to prevent BLAS-related runtime errors.
  * **Improved documentation:** Installation section now includes step-by-step configuration guide with examples for Linux and macOS, CUDA version adjustment instructions, and expanded troubleshooting section.
* **[2026/04/30]** 🧰 Extra-binner audit (`-e v` / `-e l`).
  * VAMB 4.x is pinned to Python 3.11 by bioconda, so the BASALT 3.12
    base env **cannot install it**. VAMB 5.x is required; BASALT's
    wrapper now auto-detects the installed VAMB major version at
    runtime and emits the correct CLI form (`vamb bin default …` for
    v5, legacy `vamb …` for v4 if you still have it in a separate env).
  * VAMB 5.x cannot be co-installed with the rest of BASALT's pixi env
    without breaking the PyTorch graph, so it is **not bundled**. Run
    `BASALT --check-deps -e v` to verify presence; install in a side
    env and prepend its `bin/` to `$PATH`.
  * LorBin is not on conda or PyPI yet — install it from upstream.
    Output dir / QC paths in BASALT's wrapper now use a consistent
    lowercase `_100_lorbin_genomes` (was a mix of `LorBin` and
    `lorbin`, breaking CheckM2 invocation).
  * Both wrappers now consume **sorted** BAMs (`*_DNA-*_sorted.bam`);
    they previously fed unsorted BAMs which VAMB 5.x rejects outright.
  * Removed the dead `from lib2to3.fixes import fix_buffer` import in
    `s1_autobinners.py` and `s1e_extra_binners.py` — `lib2to3` is
    deprecated in 3.12 and slated for removal in 3.13.
  * `s2`/`s6`/`s7` filename-detection chains now recognise the
    `lorbin` substring so LorBin outputs follow the same `.fa`
    rename path as MetaBAT2/VAMB.
* **[2026/04/30]** 🧹 **MetaBinner support removed.** The `-e m` flag is
  no longer accepted; use `-e v` (VAMB) and/or `-e l` (LorBin) instead,
  or omit `-e` entirely. Existing `-e m` invocations will now fail fast
  with a clear error.
* **[2026/04/30]** 🧰 Usability + rigor pass.
  * **Absolute paths** in `-a` / `-s` / `-l` / `-hf` / `-d` / `-r` / `-c`
    are now supported — BASALT auto-symlinks them into the working
    directory at startup, so you no longer need to `cd data/` first.
  * **`--workdir DIR`** controls where intermediate files live;
    **`--outdir DIR`** controls where the final output folder is
    moved on completion. Defaults preserve the historical
    cwd-everything behaviour.
  * **`--version`**, **`--check-deps`**, **`--dry-run`** flags added.
    `--check-deps` verifies that `bowtie2`, `samtools`, `minimap2`,
    `metabat2`, `checkm2` (or `checkm`), `blastn`, `spades.py`, `pilon`,
    `jgi_summarize_bam_contig_depths`, etc. are on `$PATH` before you
    submit a long job. Honours `-q` and `-e` selection.
  * **Run manifest** written to `<workdir>/BASALT_run_manifest.json`
    capturing argv, version, hostname, parameters and resolved inputs
    at start; `status`, `duration_seconds` and `error` are appended on
    completion (success *or* failure). Reproducibility-friendly.
  * **Unified logger** — every CLI/banner message now goes through
    `basalt.logger`, which dual-writes timestamped lines to stdout
    **and** `Basalt_log.txt`. **Total elapsed time** is logged at
    pipeline end, including aborted runs.
  * **Backend-aware messages** — print/log lines no longer say
    "checkm" when CheckM2 is the active backend (and vice versa).
  * **`-s` parser** now accepts absolute-path read pairs via `;` as
    the pair separator (or as a flat comma list auto-paired by twos),
    while preserving the legacy `r1.fq,r2.fq/d2_r1.fq,d2_r2.fq`
    form for backward compatibility.
  * **Bug fixes** in autobinning / refinement / reassembly format
    detection (`assembly_list[0]` → per-iter `item`), in S7 contig
    retrieval (`cutoff0`/`cutoff1` use-before-assign), in S1e
    LorBin command construction, in S9 long-read SAM cleanup, in
    S1p bin-pair logic, and several typos in user-visible strings.
* **[2026/04/28]** ⚡ Hot-path performance overhaul.
  * **Single-pass SAM parsing** in S9 reassembly, S7lr long-read
    polishing and S7p gap-filling. Each SAM is scanned once and
    per-bin FASTQ buffers are flushed in a thread pool — replaces a
    code path that opened-appended-closed each `_seq_R{1,2}.fq` file
    once per read.
  * **Coverage-matrix cache in S4** — each `Coverage_list_*.txt` is
    parsed and tokenised once instead of up to 3× per de-rep
    iteration.
  * **Parallel per-dataset mapping in S1** — bowtie2 / minimap2
    invocations across multiple samples can now run concurrently via
    the new `BASALT_DATASETS_PARALLEL` env var (default `2`). See
    [*Performance & tuning*](#-performance--tuning).
* **[2026/04/27]** 🤗 **BASALT v1.2.0** rewritten as a pip-installable Python
  package. No more `install.sh` / `chmod +x` rituals. **`pixi` is now the
  recommended installer** — a single `pixi install` builds the full
  conda + pip env and sets all activation env vars. CheckM and CheckM2
  paths unified through a runtime backend switch (`-q checkm` /
  `-q checkm2`). All step modules grouped under `basalt.steps`,
  `basalt.module`, `basalt.core`, `basalt.ml`. See *Installation* below.
* **[2025/12/16]** 🤗 BASALT v1.2.0 (preview) released under MIT.
* **[2024/06/12]** BASALT v1.1.0 released.
* **[2024/03/11]** Paper accepted at *Nature Communications*.
* **[2023/08/18]** BASALT v1.0.0 released.

## 🌉 Workflow

<img src="fig/workflow.png" style="zoom: 75%;" />

The pipeline is split into four functional **modules** (selected via
`--module`):

| Module        | What it does                                                   | Key steps |
|---------------|----------------------------------------------------------------|-----------|
| `autobinning` | Map reads, run binners (MetaBAT2 / MaxBin2 / CONCOCT / SemiBin / LorBin), within-group dereplication, multi-assembly dereplication | S1 / S1e / S1p / S2 / S3 / S4 |
| `refinement`  | DL-based outlier removal, within-group contig retrieval, optional polishing | S4 / S5 / S6 / S7 / S7lr |
| `reassembly`  | OLC elongation, short-read / hybrid reassembly, post-reassembly OLC dereplication | S4 / S8 / S9 / S9p / S10 |
| `datafeeding` | Integrate external binsets into a BASALT run, then re-run S4 + S5 | (triggered by `-d`) |

The default `--module all` runs autobinning → refinement → reassembly in
sequence, then optionally `datafeeding` if `-d` is supplied.

## 📊 Supported inputs

| Input type | Flag | Notes |
|---|---|---|
| Short-read assembly | `-a` | Files: `.fa` / `.fna` / `.fasta`, plus `.gz` / `.tar.gz` / `.zip`. Absolute or relative paths both work — BASALT auto-symlinks into the working directory at startup. |
| Long-read or hybrid assembly | `-a` | Same formats |
| Paired-end short reads | `-s` | `r1.fq,r2.fq/d2_r1.fq,d2_r2.fq` (pairs separated by `/`). For absolute paths use `;` instead of `/` between pairs (since `/` collides with directory separators), e.g. `/abs/r1_1.fq,/abs/r1_2.fq;/abs/r2_1.fq,/abs/r2_2.fq`. A flat comma list is also accepted and auto-paired by twos. |
| Long reads (ONT / PacBio CLR) | `-l` | |
| PacBio HiFi reads | `-hf` | |
| Pre-existing binsets to import | `-d` | Triggers the `datafeeding` workflow |
| Refinement on existing binset | `-r` + `-c` | Run refinement only against an existing bin folder |

> **Path note:** absolute paths are now supported in every input flag.
> BASALT will create a `<basename> -> <absolute target>` symlink in the
> working directory so the rest of the pipeline (which references files
> by basename) keeps working unchanged.

---

## ⚙️ Installation

BASALT v1.2.0 is a regular Python package. The **recommended** workflow
is [`pixi`](https://pixi.sh) — it reads the bundled `pixi.toml`, builds a
fully-pinned conda + pip environment, sets all required env vars, and
exposes BASALT as an editable install in a single step.

> ⚠️ **Python 3.12 is required.** BASALT v1.2.0 pins to Python 3.12.x —
> earlier versions are unsupported (PyTorch, scikit-learn, and several
> bioinformatics pinnings only line up cleanly on 3.12).

### ⭐ Recommended: pixi (one-shot install)

Pixi replaces steps 1–4 of the manual workflow below. You get the conda
env, the editable BASALT install, and the activation-time env vars
(`BASALT_WEIGHT`, `CHECKM2DB`) all from one declarative file.

#### Step 1: Install pixi

```bash
curl -fsSL https://pixi.sh/install.sh | sh
```

#### Step 2: Clone the repository

```bash
git clone https://github.com/EMBL-PKU/BASALT.git
cd BASALT
```

#### Step 3: Configure paths (REQUIRED before installation)

> ⚠️ **IMPORTANT:** You **MUST** edit `pixi.toml` before running `pixi install`.

Open `pixi.toml` and modify the `[activation.env]` section (around lines 78-87) to match your filesystem:

```toml
[activation.env]
BASALT_WEIGHT = "/your/actual/path/to/basalt_weights"        # CHANGE THIS
CHECKM2DB     = "/your/actual/path/to/checkm2db/CheckM2_database/uniref100.KO.1.dmnd"  # CHANGE THIS
# CHECKM_DATA_PATH = "/your/actual/path/to/checkmdb"   # CHANGE THIS if running -q checkm
```

**Example configurations:**

```toml
# Linux example
BASALT_WEIGHT = "/home/username/basalt_weights"
CHECKM2DB     = "/home/username/databases/checkm2db/CheckM2_database/uniref100.KO.1.dmnd"

# macOS example
BASALT_WEIGHT = "/Users/username/basalt_weights"
CHECKM2DB     = "/Users/username/databases/checkm2db/CheckM2_database/uniref100.KO.1.dmnd"
```

**Optional: Adjust CUDA version**

If your system uses a different CUDA version, modify the `[system-requirements]` section (around line 13):

```toml
[system-requirements]
cuda = "12"  # Change to "11" or "13" based on your CUDA version
```

Check your CUDA version with:
```bash
nvidia-smi  # Look at "CUDA Version" in the top-right
# or
nvcc --version
```

**For CPU-only systems:**

If you don't have a GPU, edit `pixi.toml`:
1. Remove or comment out the `[system-requirements]` section
2. Change `pytorch-gpu = "*"` to `pytorch-cpu = "*"` in the `[dependencies]` section

#### Step 4: Install the environment

```bash
pixi install
```

This will:
- Create a conda environment with Python 3.12
- Install all bioinformatics tools (metabat2, checkm2, bowtie2, etc.)
- Install BASALT in editable mode
- Set up environment variables

#### Step 5: Download model weights and databases

```bash
pixi run download-weights     # Downloads BASALT DL models (~100 MB) → $BASALT_WEIGHT
pixi run checkm2-db           # Downloads CheckM2 DIAMOND DB (~3 GB) → directory containing $CHECKM2DB
```

#### Use the env

Two equivalent ways to run BASALT:

```bash
# A. Persistent shell — env vars + tools active until you `exit`
pixi shell
BASALT --help

# B. One-shot — runs a single command inside the env, then exits
pixi run BASALT --help
```

`pixi shell` / `pixi run` automatically export `BASALT_WEIGHT` and
`CHECKM2DB`, so you don't need to touch `~/.bashrc`. From any working
directory you can simply call `BASALT …`.

#### Verify installation

```bash
pixi shell

# Check BASALT version
BASALT --version

# Verify all dependencies are installed
BASALT --check-deps

# Check environment variables
echo $BASALT_WEIGHT
echo $CHECKM2DB

# If using extra binners (VAMB or LorBin)
BASALT --check-deps -e v,l
```

#### Recommended directory structure

```
/your/base/directory/
├── basalt_weights/          # BASALT model weights (set in pixi.toml)
│   └── BASALT/
│       ├── model1.pth
│       └── model2.pth
├── databases/
│   └── checkm2db/
│       └── CheckM2_database/
│           └── uniref100.KO.1.dmnd  # Path set in pixi.toml
└── BASALT/                  # BASALT source code (this repository)
    ├── pixi.toml
    ├── src/
    └── ...
```

GPU note: `pixi.toml` declares `cuda = "12"` under `[system-requirements]`,
so on a CUDA 12+ host pixi will pull GPU PyTorch automatically.

---

### Alternative: manual conda + pip / uv

Use this path only if pixi isn't available on your host (e.g. some HPC
clusters block user-side installers). Steps 1–4 below reproduce what
`pixi install` does in one command.

#### 1. Create the conda environment (non-Python deps)

```bash
conda create -n basalt_env -c conda-forge -c bioconda \
    python=3.12 perl \
    metabat2 maxbin2 concoct semibin checkm2 \
    bowtie2 bwa samtools blast diamond hmmer prodigal \
    spades unicycler idba racon pilon pplacer \
    minimap2 miniasm megahit bedtools entrez-direct ncbi-vdb --yes

conda activate basalt_env
```

This pulls in `jgi_summarize_bam_contig_depths` (via metabat2),
`checkm2`, all aligners and assemblers — i.e. every external tool BASALT
shells out to.

#### 2. Install BASALT

Two interchangeable installers — pick whichever matches your workflow:

##### Option A. `pip` (always available)

```bash
# From PyPI (released)
pip install BASALT

# Or, from a clone (development / editable install)
git clone https://github.com/EMBL-PKU/BASALT.git
cd BASALT
pip install -e .
```

##### Option B. `uv` (10–100× faster than pip)

```bash
# In an activated conda env (uv reuses its Python interpreter)
uv pip install BASALT

# Or from a clone:
git clone https://github.com/EMBL-PKU/BASALT.git
cd BASALT
uv pip install -e .
```

If you don't already have [`uv`](https://github.com/astral-sh/uv):

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

`uv` resolves and installs into the *currently active* Python environment
(the conda `basalt_env` you created in step 1), so no extra setup is
needed beyond running `uv pip install …` in place of `pip install …`.

##### Option C. Pure-`uv` workflow (no conda, advanced users)

If you prefer to avoid conda entirely and supply the non-Python tools
yourself (system package manager, custom builds, etc.):

```bash
uv venv --python 3.12 basalt_env
source basalt_env/bin/activate
uv pip install BASALT
# Then ensure metabat2, checkm2, bowtie2, spades, … are on $PATH yourself.
```

---

Whichever installer you pick, the result is the same:

* Two console scripts registered on `PATH`:
  * `BASALT` — the main pipeline
  * `BASALT_models_download` — downloads the deep-learning weights
* Bundled Perl helpers (`calc.kmerfreq.pl`, `Cytoscapeviz.pl`) ship
  inside the wheel.

No `chmod`, no `install.sh`, no manual file copying.

#### 3. Download model weights

```bash
# Pick any directory you like; remember the path.
BASALT_models_download --path /path/to/basalt_weights
```

#### 4. Set environment variables

Append to `~/.bashrc` (or `~/.zshrc`):

```bash
export CHECKM2DB=/path/to/checkm2db/CheckM2_database/uniref100.KO.1.dmnd
export CHECKM_DATA_PATH=/path/to/checkmdb           # only if you use --quality-check checkm
export BASALT_WEIGHT=/path/to/basalt_weights        # the path you used in step 3
```

Then `source ~/.bashrc`.

CheckM2 and CheckM database files (and a Singularity image) are
mirrored at:

> https://drive.google.com/drive/folders/1d0e_2FpYRBAZLwKXl8fA-yDK4b5PBA_E

### Alternative: Singularity (China mainland)

Place `basalt.sif` in your home directory:

```bash
singularity run basalt.sif BASALT -a as1.fa \
    -s S1_R1.fq,S1_R2.fq/S2_R1.fq,S2_R2.fq -t 32 -m 128

# If basalt.sif is elsewhere, mount the data directory:
singularity run -B /media/emma basalt.sif BASALT -h
```

For long jobs, prefer `sbatch` on a cluster, or `screen` on a workstation:

```bash
screen -dmS my_basalt_job bash -c 'bash basalt.sh > log_basalt'
```

The image bundles `checkm`, `checkm2`, `semibin`, `bowtie2`, `bwa`, etc.;
each can be invoked directly:

```bash
singularity run basalt.sif bowtie2 -h
```

---

## 🚀 Quick start

### Minimal short-read run

```bash
BASALT \
    -a assembly.fa \
    -s r1.fq,r2.fq \
    -t 32 -m 128
```

### Multiple short-read datasets + multiple assemblies

```bash
BASALT \
    -a as1.fa,as2.fa,as3.fa \
    -s d1_r1.fq,d1_r2.fq/d2_r1.fq,d2_r2.fq/d3_r1.fq,d3_r2.fq \
    -t 60 -m 250
```

### Hybrid (short + long + HiFi)

```bash
BASALT \
    -a assembly.fa \
    -s sr1_r1.fq,sr1_r2.fq \
    -l ont1.fq,ont2.fq \
    -hf hifi1.fq \
    -t 60 -m 250
```

### Refinement of a pre-existing binset

```bash
BASALT \
    -r my_existing_binset \
    -c Coverage_matrix_for_binning_assembly.fa.txt \
    -s r1.fq,r2.fq \
    -t 32 -m 128 \
    --module refinement
```

### CheckM (legacy) instead of CheckM2

```bash
BASALT -a as.fa -s r1.fq,r2.fq -t 32 -m 128 -q checkm
```

### Custom workdir + outdir + absolute-path inputs

Run from anywhere; keep bulky intermediates separate from the final binset:

```bash
BASALT \
    -a /data/asm1.fa,/data/asm2.fa \
    -s /data/r1_1.fq,/data/r1_2.fq;/data/r2_1.fq,/data/r2_2.fq \
    -t 64 -m 128 \
    -o my_run \
    --workdir /scratch/my_basalt_work \
    --outdir  /results/my_basalt_out
```

* All process / checkpoint / log files land in `--workdir`.
* When the pipeline finishes, the `<-o>` folder is moved into `--outdir`.
* Absolute paths are accepted directly — no need to `cd` into a data directory first.
* `--workdir` defaults to the current directory; `--outdir` defaults to `--workdir` (legacy behaviour).

### Pre-flight check & dry-run

```bash
# Verify every required external tool is on $PATH (honours -q and -e):
BASALT --check-deps

# Resolve all inputs/parameters, write the run manifest, then exit:
BASALT -a asm.fa -s r1.fq,r2.fq --dry-run

# Print version:
BASALT --version
```

---

## 🎛 CLI options (summary)

Run `BASALT --help` for the authoritative list. The most-used flags:

| Flag | Meaning | Default |
|---|---|---|
| `-a` | Assembly FASTA list (comma-separated) | — |
| `-s` | Paired-end short reads | — |
| `-l` | Long reads (ONT / PacBio CLR) | — |
| `-hf` | PacBio HiFi reads | — |
| `-d` | External binsets to ingest (data feeding) | — |
| `-r` | Existing binset folder for refinement-only | — |
| `-c` | Coverage matrix file(s) | — |
| `-b` | Binset folders for standalone S4 dereplication | — |
| `-t` | Threads | `4` |
| `-m` | RAM in GB | `32` |
| `-q` / `--quality-check` | `checkm2` (default) or `checkm` | `checkm2` |
| `--module` | `all` / `autobinning` / `refinement` / `reassembly` | `all` |
| `--mode` | `continue` (resume from last checkpoint) or `new` | `continue` |
| `--sensitive` | `quick` / `sensitive` / `more-sensitive` | `sensitive` |
| `--refinepara` | `quick` / `deep` | `quick` |
| `--min-cpn` | Minimum bin completeness % | `35` |
| `--max-ctn` | Maximum bin contamination % | `20` |
| `-e` | Extra binners: `v` = VAMB (≥5; install in a side env), `l` = LorBin (install from upstream). MetaBinner removed in v1.2.x. Run `BASALT --check-deps -e v,l` to verify presence before launching. | — |
| `-o` | Output folder name | `Final_binset` |
| `-w` / `--workdir` | Directory for intermediate process files | cwd |
| `--outdir` | Directory the final output folder is moved to on completion | same as `--workdir` |
| `-V` / `--version` | Print BASALT version and exit | — |
| `--check-deps` | Verify required external tools are on `$PATH` and exit (non-zero if any missing) | — |
| `--dry-run` | Resolve config + write run manifest, then exit before dispatch | — |

---

## 📑 Logging & reproducibility

Every BASALT run produces three machine-readable artefacts in `--workdir`:

| File | What it contains |
|---|---|
| `Basalt_log.txt` | Timestamped pipeline log (also streamed to stdout). Driven by the unified `basalt.logger` module — every message is prefixed with `[YYYY-MM-DD HH:MM:SS] [LEVEL]`. |
| `Basalt_checkpoint.txt` | Resume points used by `--mode continue`. |
| `BASALT_run_manifest.json` | Full reproducibility manifest. Captures BASALT version, complete `argv`, parsed parameters, normalised inputs, hostname / user / pid / Python version, start time, and on completion: `status` (`success` / `aborted`), `duration_seconds`, `ended_at`, and `error` (if applicable). |

The total elapsed wall-clock is logged at end of every run (success
*and* failure) as `BASALT pipeline finished in H:MM:SS (123.4 s)` —
useful for cluster accounting and benchmarking.

---

## ⚡ Performance & tuning

Most of BASALT's tuning happens via CLI flags (`-t`, `-m`,
`--sensitive`, `--refinepara`). One additional knob lives in an
**environment variable**:

| Env var | Default | What it controls |
|---|---|---|
| `BASALT_DATASETS_PARALLEL` | `2` | Number of per-dataset bowtie2 / minimap2 mapping jobs to run concurrently inside S1's mapping helpers (`mapping`, `mapping_lr_o`, `mapping_hifi_split`, `mapping_hifi_minimap`). Each worker is given `-t // BASALT_DATASETS_PARALLEL` mapper threads, so total CPU stays roughly bounded. RAM scales linearly (each mapper holds its own index in memory). |

The clamp is `workers = max(1, min(env_var, num_samples))`, so
setting it higher than your sample count is harmless — it just runs
fewer workers.

### When to scale `BASALT_DATASETS_PARALLEL`

| Number of input samples | Recommended | Effect |
|---|---|---|
| 1 sample (any data type) | any value | Auto-clamps to 1 — no concurrency. |
| 2 samples | `2` (default) | 2 mappers in parallel, half threads each. |
| 3-4 samples | `3` | |
| 5+ samples | `4` | Beyond 4 the win disappears: bowtie2 / minimap2 thread efficiency drops fast below ~8 threads each. |

Example for a 4-sample run on a 64-thread node with plenty of RAM:

```bash
BASALT_DATASETS_PARALLEL=4 \
  BASALT -a as.fa -s s1_R1.fq,s1_R2.fq/s2_R1.fq,s2_R2.fq/s3_R1.fq,s3_R2.fq/s4_R1.fq,s4_R2.fq \
         -t 64 -m 256 -o my_run
```

### Notable internal speedups (no flags required)

These changes are automatic — no tuning needed, listed for reference:

* **Single-pass SAM parsing** — `parse_sam`, `parse_sam_bwa` and
  `parse_lr_sam` in `basalt/steps/s9_reassembly.py`,
  `s7lr_finding_sr_contigs.py` and `s7p_gap_filling.py` now do **one
  scan** over each SAM file and flush per-bin FASTQ buffers via
  `concurrent.futures.ThreadPoolExecutor`. Replaces an old code path
  that opened-appended-closed each output FASTQ file once per
  alignment record.
* **Coverage-matrix cache** — `basalt/steps/s4_multiple_assembly_comparator.py`
  caches each `Coverage_list_*.txt` parse via a process-lifetime
  `_parse_coverage(path)` helper, so callers
  (`genome_contigs_recorder`, `binset_comparator`,
  `new_selected_bins_generator`, `record_bin_coverage`) reuse a
  single tokenised view instead of re-reading the file from disk.

Both are designed to be output-preserving. Verified equivalent on
fixture data; if you're migrating from a pre-2026/04/28 BASALT and
want to be sure, run the old and new versions on the same input and
compare `OLC_quality_report.tsv` side-by-side.

---

## 📦 Package architecture (for developers)

```
basalt/
├── cli.py                 # `BASALT` entry point + dependency check + run manifest
├── logger.py              # unified stdout + Basalt_log.txt logging facade
├── qc_backend.py          # CheckM2 / CheckM unified abstraction
├── shell.py               # subprocess wrappers
├── core/                  # cross-step utilities (cleanup / data_feeding / final_drep)
├── module/                # 4 phase entry points — names match `--module`
│   ├── autobinning.py
│   ├── refinement.py
│   ├── reassembly.py
│   └── datafeeding.py
├── steps/                 # 17 numbered pipeline steps S1–S10
├── ml/                    # deep-learning model + weights downloader
└── scripts/               # bundled Perl helpers (calc.kmerfreq.pl etc.)
```

`basalt/module/<phase>.py` configures the QC backend at function entry
based on `QC_software`, then dispatches to `basalt/steps/*` modules.
A single set of step modules serves both CheckM2 and CheckM via the
`qc_backend` abstraction.

See [`src/basalt/ARCHITECTURE.md`](src/basalt/ARCHITECTURE.md) for an
in-depth tour.

---

## 🛠 Troubleshooting

### `BASALT: command not found`

* **pixi users:** make sure you're inside `pixi shell` (or wrapping the
  call with `pixi run`). `which BASALT` should point into
  `<repo>/.pixi/envs/default/bin/BASALT`.
* **conda + pip users:** make sure you activated the conda environment
  **and** ran `pip install BASALT` (or `pip install -e .` from source)
  inside it. `which BASALT` should point into your conda env.

```bash
which BASALT          # should point into the active env
pip show BASALT       # should print v1.2.0
```

### `ERROR: DIAMOND database not found` (CheckM2)

**For pixi users:** Make sure you edited `pixi.toml` to set the correct `CHECKM2DB` path before running `pixi install`. Then run:

```bash
pixi run checkm2-db
```

**For manual installation:** Download the database:

```bash
checkm2 database --download --path /path/to/checkm2db
export CHECKM2DB=/path/to/checkm2db/CheckM2_database/uniref100.KO.1.dmnd
```

### Environment variables not set

**For pixi users:** If `echo $BASALT_WEIGHT` or `echo $CHECKM2DB` returns empty:

1. Check that you edited `pixi.toml` before running `pixi install`
2. Make sure you're inside `pixi shell` or using `pixi run`
3. Re-run `pixi install` after editing `pixi.toml`

**For manual installation:** Add to `~/.bashrc` or `~/.zshrc`:

```bash
export BASALT_WEIGHT=/path/to/basalt_weights
export CHECKM2DB=/path/to/checkm2db/CheckM2_database/uniref100.KO.1.dmnd
```

Then `source ~/.bashrc` (or `~/.zshrc`).

### CUDA version mismatch

If you see PyTorch CUDA errors:

1. Check your CUDA version: `nvidia-smi` or `nvcc --version`
2. Edit `pixi.toml` line 13 to match your CUDA version:
   ```toml
   cuda = "11"  # or "12" or "13"
   ```
3. Re-run `pixi install`

### `libopenblas` or BLAS-related errors

The updated `pixi.toml` includes `libopenblas` as a dependency. If you still encounter BLAS errors:

```bash
# For pixi users
pixi install  # Re-run to ensure libopenblas is installed

# For manual conda users
conda install -c conda-forge libopenblas
```

### `samtools: error while loading shared libraries: libcrypto.so.1.0.0`

Symlink the newer libcrypto:

```bash
cd "$CONDA_PREFIX/lib"
ln -s libcrypto.so.1.1 libcrypto.so.1.0.0
```

### Missing external tools

If a run dies hours in with errors like `bowtie2: command not found`
or `checkm2 ... cannot find database`, run the pre-flight check
**before** submitting next time:

```bash
BASALT --check-deps
```

It prints a green `[ OK ]` / red `[MISS]` line for every executable
BASALT will shell out to (honours `-q checkm` vs `checkm2` and
`-e m`/`v`/`l` for extra binners).

### `FileNotFoundError: 'quality_report.tsv'` in `Contig_recruiter_main`

The current step produced too few bins for CheckM2 to score, usually
because read coverage is very low for the dataset. Increase coverage,
combine with another sample, or relax `--min-cpn` / `--max-ctn`.

### Long-read–only mode fails

LRS-only mode (Nanopore / PacBio CLR without short reads) is
**unsupported in v1.2.0** — supply at least one short-read pair via `-s`,
or use HiFi-only via `-hf`.

---

## 🤝 FAQ

> **Q. Same contigs are used in both single-assembly and co-assembly modes — does this affect output?**
> A. Redundant bins from overlapping assemblies are detected and removed
> in S3 (within-group) and S4 (cross-assembly) dereplication. Final
> output is non-redundant.

> **Q. BASALT takes longer than metaWRAP — can I speed it up?**
> A. Yes:
> ```bash
> --sensitive quick --refinepara quick --module autobinning
> ```
> trades some MAG quality / count for shorter runtime. BASALT is also
> single-pass over multiple assemblies, so on a multi-assembly project
> it is often *faster* than running metaWRAP per assembly + dRep.
>
> If you have multiple input samples and RAM headroom, also set
> `BASALT_DATASETS_PARALLEL=3` (or `4`) to run several
> bowtie2 / minimap2 jobs concurrently in S1 — see
> [*Performance & tuning*](#-performance--tuning).

> **Q. How do I refine an existing bin set?**
> A. Use `--module refinement` with `-r my_binset -c coverage_matrix.txt`.
> Or import bins via `-d my_binset_folder` for the data-feeding workflow.

> **Q. Where do results go?**
> A. By default, into a `<-o>` folder under the current working directory.
> Use `--workdir DIR` to direct intermediate files (logs, checkpoints,
> per-step scratch) to a chosen folder, and `--outdir DIR` to have the
> final binset moved to a separate location on completion. The legacy
> "everything in cwd" layout is preserved when neither flag is set.

> **Q. How do I reproduce a previous run?**
> A. Look at `BASALT_run_manifest.json` in the workdir — it captures
> the exact `argv`, BASALT version, parameters, hostname and timing of
> the run. The same command line can be replayed; combine with
> `--mode new` to start a fresh project.

---

## 📜 License

MIT — see [LICENSE](LICENSE) for the full text.

## ✏️ Citation

If you use BASALT in published research, please cite:

> Z. Qiu, L. Yuan, C. Lian, B. Lin, J. Chen, R. Mu, X. Qiao, L. Zhang,
> Z. Xu, L. Fan, Y. Zhang, S. Wang, J. Li, H. Cao, B. Li, B. Chen,
> C. Song, Y. Liu, L. Shi, Y. Tian, J. Ni, T. Zhang, J. Zhou, W. Zhuang,
> K. Yu. **BASALT refines binning from metagenomic data and increases
> resolution of genome-resolved metagenomic analysis.**
> *Nature Communications* **15**, 2179 (2024).
> https://doi.org/10.1038/s41467-024-46539-7

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

## 📚 References

1. Qiu, Z. et al. *Nature Communications* **15**, 2179 (2024). [BASALT]
2. Sieber, C. M. K. et al. *Nature Microbiology* **3**, 836–843 (2018). [DAS_Tool]
3. Uritskiy, G. V. et al. *Microbiome* **6**, 158 (2018). [metaWRAP]
4. Olm, M. R. et al. *The ISME Journal* **11**, 2864–2868 (2017). [dRep]
5. Xue, W. et al. *Nature Communications* **16**, 9353 (2025). [LorBin]

---

For bug reports, feature requests, or compilation issues, please open
an issue on GitHub or contact <yuke.sz@pku.edu.cn>.
