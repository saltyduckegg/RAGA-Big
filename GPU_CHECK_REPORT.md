# GPU Acceleration Verification Report

## Summary
The software used in the existing workflow was examined to determine if GPU acceleration is enabled. The analysis found that while **racon** supports GPU acceleration (via CUDA), it is **not** enabled in the current configuration and scripts.

## Findings

### 1. Software with Potential GPU Support
Among the tools listed in the `workflow.txt` and used in the `bin/` scripts, **racon** is the only tool that natively supports GPU acceleration for its core tasks (polishing).

Other tools such as `minimap2`, `hifiasm`, `nucmer` (MUMmer), `RagTag`, `bedtools`, `samtools`, `bwa`, `pilon`, and `spades.py` are primarily CPU-based in the versions and configurations used here.

### 2. Script Analysis
The following scripts invoke `racon`:
- `bin/RAGA-same.sh`
- `bin/RAGA-diff.sh`

In all instances, `racon` is called with the following pattern:
```bash
racon <reads> <overlaps> <target> -t $thr > <output>
```
The `-t` flag specifies the number of threads (CPU). There are **no** flags enabling GPU acceleration, such as:
- `-g` (or similar, depending on version)
- `-c <batches>` (used to control CUDA batches in Racon)

### 3. Installation and Environment
- The `README.md` installation instructions for `racon` use standard `cmake` without the `-Dracon_enable_cuda_c=ON` flag, which is required to build `racon` with CUDA support.
- No environment variables (e.g., `CUDA_VISIBLE_DEVICES`) are set in the scripts to facilitate GPU usage.

## Conclusion
**GPU acceleration is NOT enabled** in the existing process. To enable it, one would need to:
1. Recompile `racon` with CUDA support (`-Dracon_enable_cuda_c=ON`).
2. Update the scripts (`bin/RAGA-same.sh` and `bin/RAGA-diff.sh`) to include GPU-specific flags (if required by the specific version/build of racon, though typically the CUDA build might use it by default or requires `-c`).
