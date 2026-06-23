# How to run the AF3 pairwise co-folding pipeline

This pipeline takes a FASTA file with N proteins, predicts the 3D structure of **every possible pair** (every heterodimer plus each protein paired with itself = homodimer), and gives you confidence metrics for every prediction. It uses AlphaFold 3 under the hood and runs on the EMBL HPC cluster.

For N proteins, you get N(N+1)/2 predictions. So 10 proteins → 55 pairs, 12 proteins → 78 pairs, etc.

---

## Before you start

You need to be logged into an EMBL **login node** (`login1.cluster.embl.de` or similar). This is just a normal Linux machine where you submit jobs to the cluster — you don't run the actual compute here.

The cluster uses a job scheduler called **SLURM**. You write a job (or a script does it for you), submit it, and SLURM finds a free CPU/GPU somewhere in the cluster, runs your job there, and saves the output. Your jobs run independently of your terminal — **you can close your laptop and come back the next day**.

You **don't** need to set up python or load any modules yourself — the launcher scripts (`run_pipeline.sh`, `multifold.sh`, `run_qc.sh`) auto-load the right AlphaFold3 module on first call. If you've already loaded something else, that's fine, the scripts only top up what's missing.

---

## Run it

One command:

```bash
/g/typas/Personal_Folders/Nic/sophie_viral_cofolding/scripts/run_pipeline.sh \
    /path/to/your_proteins.faa \
    /g/typas/Personal_Folders/Nic/sophie_viral_cofolding/runs/<your_run_name>
```

Convention: name run dirs `YYYY-MM-DD_<short_description>` under the project's `runs/` folder (e.g. `runs/2026-06-23_dimer_all_pairs/`). Past runs live there for reference — see git log for what went in each.

The script will:
1. Create the output directory
2. Copy your FASTA file and the pipeline scripts into it (so you have a frozen snapshot of exactly what produced your results)
3. Generate one input JSON per pair
4. Submit two SLURM jobs per pair: one for the MSA, one for the structure prediction
5. Print a status command and a QC command for you to use later

Then it exits. Your jobs are now in SLURM's queue.

**Why two jobs per pair?** AF3 has two stages. The first (called "MSA") searches huge sequence databases for related proteins — this is CPU work, not GPU. The second is the actual structure prediction — this needs a GPU. Running them as one big GPU job would waste the GPU during the CPU phase. The two are chained automatically: the GPU job only starts after the CPU job finishes successfully.

---

## How long does it take?

Per pair, roughly:
- MSA (CPU): 10–15 min
- Structure prediction (GPU): 10–60 min, longer for bigger pairs

Many pairs run **in parallel**, so total wall time is usually 1–4 hours regardless of how many pairs. The main bottleneck is GPU availability — sometimes you wait in queue.

---

## Watch progress

The pipeline prints these commands when it submits. You can run them at any time:

**How many of your jobs are left in the queue?**
```bash
squeue -u $USER -h | awk '$3 ~ /^msa_|^inf_/' | wc -l
```
When this prints `0`, everything is done.

**Per-job status:**
```bash
squeue -u $USER | grep -E "msa_|inf_"
```
The `ST` column means:
- `R` = currently running
- `PD` = pending (waiting in queue, this is normal)
- (job disappears from the list once it finishes)

**Did anything fail?**
```bash
cat /path/to/where/you/want/output/logs/jobid_map.tsv
# Take any jobid from there:
sacct -j <jobid> --format=JobID,JobName,State,Elapsed,ExitCode
```
`State=COMPLETED` is what you want. `FAILED`, `TIMEOUT`, `OUT_OF_MEMORY` mean look at the log file in `<output>/logs/`.

---

## Once everything is done: run QC

The pipeline tells you the exact command. It looks like:

```bash
/path/to/where/you/want/output/scripts/run_qc.sh /path/to/where/you/want/output
```

This wrapper loads the AlphaFold3 module (for python + pandas/numpy/matplotlib) and runs `03_qc.py` against your output directory.

This produces, in `<output>/qc/`:

| file | what it is |
|---|---|
| `qc_summary.tsv` | one row per pair with all metrics. Open in Excel or pandas. |
| `iptm_heatmap.pdf` | **the headline figure**: N×N grid showing interface confidence for every pair |
| `ranking_score_heatmap.pdf` | same grid but with AF3's combined score |
| `pae_<pair>.png` | one error plot per pair (lower = AF3 more confident about position) |
| `interface_contacts.tsv` | number of inter-chain contacts and closest approach distance per pair |

---

## Where are the actual structures?

For each pair, in `<output>/results/<pair>_<timestamp>/`:
- `<pair>_model.cif` — the predicted structure. Open in PyMOL or ChimeraX.
- 5 subfolders `seed-1_sample-{0..4}/` with all 5 alternative samples AF3 generated

You'll usually only look at `<pair>_model.cif` (the best-ranked sample). The others are useful if the best one has a clash flag and you want to see if a different sample is clean.

---

## Reading the confidence metrics (the 10-second guide)

Look at these in this order:

1. **`ipTM`** (0–1): how confident AF3 is about how the chains sit relative to each other.
   - ≥ 0.8 — strong, real interaction
   - 0.6 – 0.8 — medium
   - < 0.5 — likely noise, no real interaction predicted

2. **`has_clash`** (0 or 1): if 1, the structure has overlapping atoms (physically impossible). **Ignore the ipTM in that case** — the prediction is broken.

3. **`ranking_score`**: AF3's combined score, used to pick the best of 5 samples.
   - ~0.9 — excellent
   - around 0 — bad
   - **anywhere near −99** — clash penalty fired, see `has_clash`

4. **`fraction_disordered`** (0–1): how floppy the predicted structure is. If > 0.5, the "interaction" might just be two flexible tails wagging at each other rather than a real binding interface.

5. **`iface_n_contacts`**: how many residues are in contact (<8 Å between chains). High contacts + high ipTM + no clash = strong predicted interface.

---

## Customisation

You usually don't need to change anything. But if you want to:

| You want to… | Set this environment variable before running |
|---|---|
| Force every pair onto the big H100 GPU (slower queue, more memory) | `TOKEN_H100_THRESHOLD=0` |
| Force every pair onto the regular A100 GPU | `TOKEN_H100_THRESHOLD=99999` |

Example:
```bash
TOKEN_H100_THRESHOLD=0 /g/typas/.../run_pipeline.sh proteins.faa /my/output
```

By default, pairs with a combined size > 2500 tokens (≈ residues) automatically get routed to the bigger GPU.

---

## FAQ

**Can I close my terminal?**
Yes. Once `run_pipeline.sh` exits, your jobs are in SLURM and run independently. Come back later and run the status command.

**My jobs are stuck in PD (pending). Is something broken?**
Probably not — the cluster is just busy. PD jobs are waiting their turn. If a job is PD for many hours during quiet periods, run `squeue -u $USER --start` to see when SLURM estimates it'll start.

**A job failed — what do I do?**
Look at the log: `<output>/logs/<jobname>_<jobid>.out`. The most common failure is `OUT_OF_MEMORY` on a GPU job, which means the pair is too big for the assigned GPU. Either re-run that one pair with `TOKEN_H100_THRESHOLD=0`, or skip it.

**Can I add more proteins to an existing run?**
Cleanest is to start a new output directory. The pair generation is all-vs-all, so adding even one protein means many new pairs that need to be predicted.

**What's a "token"?**
For plain proteins, one token ≈ one amino acid. AF3 talks in tokens because it can also handle small molecules and nucleic acids; for our use, total tokens of a pair = sum of the two protein lengths.

**Where is the pipeline code?**
[`/g/typas/Personal_Folders/Nic/sophie_viral_cofolding/scripts/`](scripts/) — every output directory also gets its own snapshot of this code in `<output>/scripts/`, so you can reproduce a result years later.

**Where are past runs?**
[`runs/`](runs/) holds previous co-folding jobs, named `YYYY-MM-DD_<description>`. Each has its own snapshot of inputs (`*.faa`), code (`scripts/`), AF3 outputs (`results/`), and QC (`qc/`). The git log explains what each run was for.

## Multi-chain co-fold (a single complex of N chains, not all pairs)

If you want to predict ONE specific complex (e.g. Ham1×2 + NIa-Pro as a trimer) rather than every-vs-every pairs, use `multifold.sh` instead. It takes the same FASTA / ROOT arguments:

```bash
/g/typas/.../scripts/multifold.sh <fasta_with_N_chains.faa> <runs/your_run_name>
```

For homodimers/homotrimers/etc, just put the same protein in the FASTA the appropriate number of times. The script auto-detects identical sequences and tells AF3 to share their MSA (faster, and gives AF3 better paired-MSA signal).
