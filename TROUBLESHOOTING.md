# BASALT-Air Troubleshooting Guide

## Common Issues and Solutions

### 1. GLIBC Version Error with CheckM2

**Error message:**
```
ImportError: /lib64/libpthread.so.0: version GLIBC_2.18' not found 
(required by .../libnccl.so.2)
```

**Cause:** Your system's GLIBC version is older than what CheckM2's TensorFlow dependency requires (GLIBC 2.18+).

**Solution:** The `sysroot_linux-64` package is now included in `pixi.toml` (v1.0.0+) to provide compatible GLIBC runtime libraries. To apply the fix:

```bash
# Remove existing environment
pixi clean

# Reinstall with updated dependencies
pixi install

# Verify the fix
pixi run check-deps
```

**Alternative solutions if the above doesn't work:**

1. **Use Singularity/Docker container** (recommended for HPC clusters):
   ```bash
   # Build container with all dependencies
   singularity build basalt.sif docker://condaforge/mambaforge
   ```

2. **Check your system GLIBC version:**
   ```bash
   ldd --version
   ```
   If it shows < 2.17, contact your system administrator about upgrading.

3. **Use CheckM instead of CheckM2** (temporary workaround):
   ```bash
   BASALT -a assembly.fa -s reads.fq -q checkm ...
   ```
   Note: CheckM is slower and less accurate than CheckM2.

---

### 2. CUDA Version Mismatch

**Error message:**
```
RuntimeError: CUDA error: no kernel image is available for execution on the device
```

**Solution:** Edit `pixi.toml` line 13 to match your CUDA version:

```toml
[system-requirements]
cuda = "11"  # Change to "11", "12", or "13" based on your system
```

Check your CUDA version:
```bash
nvidia-smi  # Look at "CUDA Version" in the top right
```

Then reinstall:
```bash
pixi clean
pixi install
```

---

### 3. Out of Memory (OOM) Errors

**Symptoms:**
- Process killed during binning or assembly
- "Killed" message in logs

**Solutions:**

1. **Reduce thread count:**
   ```bash
   BASALT -a assembly.fa -s reads.fq -t 8 -m 64 ...  # Lower -t and -m
   ```

2. **Use fewer binners:**
   ```bash
   BASALT -a assembly.fa -s reads.fq -b metabat2,maxbin2 ...  # Skip CONCOCT/SemiBin
   ```

3. **Process samples separately:**
   ```bash
   # Instead of: -s sample1.fq,sample2.fq,sample3.fq
   # Run individually:
   BASALT -a assembly.fa -s sample1.fq -o run1 ...
   BASALT -a assembly.fa -s sample2.fq -o run2 ...
   ```

---

### 4. Missing Database Paths

**Error message:**
```
CheckM2 database not found at: /path/to/checkm2db/...
```

**Solution:** Edit `pixi.toml` lines 85-87 with your actual paths:

```toml
[activation.env]
BASALT_WEIGHT = "/your/actual/path/to/basalt_weights"
CHECKM2DB     = "/your/actual/path/to/checkm2db/CheckM2_database/uniref100.KO.1.dmnd"
```

After editing, reload the environment:
```bash
pixi shell-hook | source
```

Or restart your shell session.

---

### 5. Permission Denied Errors

**Error message:**
```
PermissionError: [Errno 13] Permission denied: '/path/to/output'
```

**Solutions:**

1. **Check directory permissions:**
   ```bash
   ls -ld /path/to/output
   chmod 755 /path/to/output  # If you own it
   ```

2. **Use a different output directory:**
   ```bash
   BASALT -a assembly.fa -s reads.fq --outdir ~/basalt_output ...
   ```

3. **On shared clusters, use your home or scratch space:**
   ```bash
   BASALT -a assembly.fa -s reads.fq \
     --workdir $HOME/basalt_work \
     --outdir $HOME/basalt_output ...
   ```

---

### 6. Dependency Check Failures

**Command:**
```bash
pixi run check-deps
```

**Common issues:**

| Tool Missing | Solution |
|--------------|----------|
| VAMB | Install separately: `conda create -n vamb_env -c bioconda 'vamb>=5'` |
| LorBin | Clone from GitHub: `git clone https://github.com/anuradhawick/LorBin` |
| Legacy CheckM | Download database and set `CHECKM_DATA_PATH` in `pixi.toml` |

---

### 7. Slow Performance

**Optimization tips:**

1. **Use more threads** (if you have the cores):
   ```bash
   BASALT -a assembly.fa -s reads.fq -t 32 ...
   ```

2. **Skip reassembly** if you only need bins:
   ```bash
   BASALT -a assembly.fa -s reads.fq --module autobinning,refinement ...
   ```

3. **Use faster binners only:**
   ```bash
   BASALT -a assembly.fa -s reads.fq -b metabat2,semibin ...  # Skip MaxBin2/CONCOCT
   ```

4. **Pre-filter short contigs:**
   ```bash
   # Filter contigs < 1500 bp before running BASALT
   seqkit seq -m 1500 assembly.fa > assembly_filtered.fa
   ```

---

### 8. Pixi Installation Issues

**Error:** `pixi: command not found`

**Solution:**
```bash
# Reinstall pixi
curl -fsSL https://pixi.sh/install.sh | sh

# Add to PATH (add to ~/.bashrc or ~/.zshrc)
export PATH="$HOME/.pixi/bin:$PATH"

# Reload shell
source ~/.bashrc  # or source ~/.zshrc
```

---

## Getting Help

If your issue isn't covered here:

1. **Check logs:** Look in `<workdir>/logs/` for detailed error messages
2. **Run with verbose mode:** Add `-v` flag for more output
3. **Report issues:** [GitHub Issues](https://github.com/EMBL-PKU/BASALT-Air/issues)

When reporting issues, please include:
- BASALT version: `pixi run version`
- System info: `uname -a` and `ldd --version`
- Full error message and command used
- Relevant log files from `<workdir>/logs/`
