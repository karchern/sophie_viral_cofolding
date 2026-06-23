#!/usr/bin/env python3
"""QC pipeline for AF3 pair predictions.

Outputs (in qc/):
  - qc_summary.tsv             per-pair status + confidence metrics
  - iptm_heatmap.pdf           N x N ipTM heatmap
  - ranking_score_heatmap.pdf  N x N ranking_score heatmap
  - iptm_vs_ranking_scatter.pdf  ipTM vs ranking_score scatter, outliers labeled
  - pae_<pair>.png             one PAE map per successful pair
  - interface_contacts.tsv     inter-chain CA-CA contacts <8 A, min distance

AF3 ranking_score formula (verified against alphafold3/model/confidences.py
`get_ranking_score()`, github.com/google-deepmind/alphafold3):

    ranking_score = 0.8 * iptm
                  + 0.2 * ptm
                  + 0.5 * fraction_disordered
                  - 100 * has_clash

Constants in the AF3 source: _IPTM_WEIGHT=0.8, _FRACTION_DISORDERED_WEIGHT=0.5,
_CLASH_PENALIZATION_WEIGHT=100.0. Note that the disorder term is ADDED (not a
penalty): AF3 boosts predictions with high fraction_disordered so it does not
penalize itself for honestly reporting disorder.
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(os.environ.get(
    "VCO_ROOT",
    "/g/typas/Personal_Folders/Nic/sophie_viral_cofolding",
))
RESULTS = ROOT / "results"
QCDIR = ROOT / "qc"
MANIFEST = ROOT / "data" / "pair_manifest.tsv"


def find_af3_output(pair_dir: Path) -> Path | None:
    """AF3 inference writes a sibling dir named <pair>_<timestamp>/ that holds
    the *_summary_confidences.json + model.cif. Find the most recent one."""
    pair = pair_dir.name
    siblings = sorted(pair_dir.parent.glob(f"{pair}_2*"), reverse=True)
    for s in siblings:
        if (s / f"{pair}_summary_confidences.json").exists():
            return s
    # Some layouts write directly into the pair_dir itself
    if (pair_dir / f"{pair}_summary_confidences.json").exists():
        return pair_dir
    return None


def load_summary(out_dir: Path, pair: str) -> dict | None:
    p = out_dir / f"{pair}_summary_confidences.json"
    return json.loads(p.read_text()) if p.exists() else None


def load_full_confidences(out_dir: Path, pair: str) -> dict | None:
    p = out_dir / f"{pair}_confidences.json"
    return json.loads(p.read_text()) if p.exists() else None


def find_model_cif(out_dir: Path, pair: str) -> Path | None:
    for cand in [out_dir / f"{pair}_model.cif", out_dir / "model.cif"]:
        if cand.exists():
            return cand
    cifs = list(out_dir.glob("*.cif"))
    return cifs[0] if cifs else None


def parse_cif_ca_coords(cif_path: Path):
    """Return dict: chain_id -> list[(res_idx, np.array([x,y,z]))]. Minimal CIF parser
    that handles AF3's mmCIF atom_site loop without needing biopython."""
    in_loop = False
    headers: list[str] = []
    rows: list[list[str]] = []
    with open(cif_path) as f:
        for line in f:
            line_s = line.rstrip()
            if line_s == "loop_":
                in_loop = True
                headers = []
                rows = []
                continue
            if in_loop:
                if line_s.startswith("_atom_site."):
                    headers.append(line_s.split(".", 1)[1])
                    continue
                if headers and line_s and not line_s.startswith("_") and not line_s.startswith("#"):
                    rows.append(line_s.split())
                    continue
                if headers and (line_s.startswith("#") or line_s.startswith("loop_")):
                    break  # finished atom_site loop
    if not headers:
        return {}

    try:
        i_atom_id = headers.index("label_atom_id")
        i_chain   = headers.index("label_asym_id")
        i_seq     = headers.index("label_seq_id")
        i_x       = headers.index("Cartn_x")
        i_y       = headers.index("Cartn_y")
        i_z       = headers.index("Cartn_z")
    except ValueError:
        return {}

    chains: dict[str, list[tuple[int, np.ndarray]]] = {}
    for row in rows:
        if len(row) <= max(i_atom_id, i_chain, i_seq, i_x, i_y, i_z):
            continue
        if row[i_atom_id].strip('"') != "CA":
            continue
        try:
            seq = int(row[i_seq])
            xyz = np.array([float(row[i_x]), float(row[i_y]), float(row[i_z])])
        except ValueError:
            continue
        chains.setdefault(row[i_chain], []).append((seq, xyz))
    return chains


def interface_metrics(cif_path: Path, cutoff: float = 8.0) -> dict:
    chains = parse_cif_ca_coords(cif_path)
    if len(chains) < 2:
        return {"n_chains": len(chains), "n_contacts": np.nan, "min_dist": np.nan}
    chain_ids = sorted(chains.keys())
    c1 = np.stack([xyz for _, xyz in chains[chain_ids[0]]])
    c2 = np.stack([xyz for _, xyz in chains[chain_ids[1]]])
    # NxM pairwise distances
    d = np.linalg.norm(c1[:, None, :] - c2[None, :, :], axis=-1)
    return {
        "n_chains": len(chains),
        "n_contacts": int((d < cutoff).sum()),
        "min_dist": float(d.min()),
    }


def plot_metric_scatter(df: pd.DataFrame, out_pdf: Path) -> None:
    """Scatter ipTM vs ranking_score, label outliers.

    Outlier definition: residual from a robust linear fit (over clash-free points)
    exceeding 2 * MAD. Plus any clashed points are always labelled since their
    ranking_score is crushed by ~100 by the AF3 clash penalty.

    Reference line plotted: ranking_score = 0.8 * ipTM (the ipTM-only contribution,
    per AF3's get_ranking_score formula). Real points sit ABOVE this line by
    approximately 0.2 * ptm + 0.5 * fraction_disordered (when not clashing).
    """
    sub = df.dropna(subset=["iptm", "ranking_score"]).copy()
    if len(sub) < 3:
        return
    clashed = sub[sub["has_clash"] == 1.0]
    clean = sub[sub["has_clash"] != 1.0]

    fig, ax = plt.subplots(figsize=(8, 6.5))

    # Reference: pure ipTM contribution (intercept = 0, slope = 0.8)
    xs = np.linspace(0, 1, 100)
    ax.plot(xs, 0.8 * xs, "--", color="gray", alpha=0.6,
            label="ranking_score = 0.8·ipTM (no pTM, no disorder term)")

    # Robust-ish fit over clean points (just OLS; few enough points that
    # MAD-based outlier detection on residuals is fine).
    fit_text = ""
    if len(clean) >= 2:
        slope, intercept = np.polyfit(clean["iptm"], clean["ranking_score"], 1)
        ax.plot(xs, slope * xs + intercept, "-", color="steelblue", alpha=0.8,
                label=f"linear fit (clean pts): y = {slope:.2f}·x + {intercept:+.2f}")
        residuals = clean["ranking_score"].values - (slope * clean["iptm"].values + intercept)
        mad = np.median(np.abs(residuals - np.median(residuals)))
        thresh = max(0.05, 2 * 1.4826 * mad)  # 1.4826 makes MAD a std proxy
        outliers = clean[np.abs(residuals) > thresh]
        fit_text = (f"linear fit y = {slope:.3f}·x + {intercept:+.3f}\n"
                    f"outliers labelled: |residual| > {thresh:.3f} (≈2·MAD)")
    else:
        outliers = clean.iloc[0:0]

    # Plot clean points, colored by fraction_disordered
    sc = ax.scatter(clean["iptm"], clean["ranking_score"],
                    c=clean["fraction_disordered"] if "fraction_disordered" in clean else "C0",
                    cmap="viridis", vmin=0, vmax=1, s=55, edgecolor="black", linewidth=0.4,
                    label="clean (no clash)")
    if "fraction_disordered" in clean:
        cb = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
        cb.set_label("fraction_disordered")

    # Plot clashed points distinctly + always label
    if len(clashed) > 0:
        ax.scatter(clashed["iptm"], clashed["ranking_score"],
                   marker="X", s=120, color="red", edgecolor="black", linewidth=0.8,
                   label="has_clash=1 (−100 penalty)")
        for _, r in clashed.iterrows():
            ax.annotate(f"  {r['pair']}", (r["iptm"], r["ranking_score"]),
                        fontsize=8, color="red", weight="bold",
                        ha="left", va="center")

    # Label residual outliers
    for _, r in outliers.iterrows():
        ax.annotate(f"  {r['pair']}", (r["iptm"], r["ranking_score"]),
                    fontsize=7, color="black", ha="left", va="center")

    ax.set_xlabel("ipTM")
    ax.set_ylabel("ranking_score")
    ax.set_title("ipTM vs ranking_score, outliers labelled")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=8)

    # Caption / formula box (bottom of figure)
    caption = (
        "ranking_score = 0.8·ipTM + 0.2·pTM + 0.5·fraction_disordered − 100·has_clash\n"
        "(AF3 get_ranking_score, github.com/google-deepmind/alphafold3 "
        "src/alphafold3/model/confidences.py)"
    )
    if fit_text:
        caption = fit_text + "\n" + caption
    fig.text(0.5, -0.02, caption, ha="center", va="top", fontsize=7.5,
             family="monospace", wrap=True)

    plt.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def plot_pae(pae: np.ndarray, pair: str, out_png: Path, chain_breaks=None):
    fig, ax = plt.subplots(figsize=(5, 4.5))
    im = ax.imshow(pae, cmap="bwr", vmin=0, vmax=30, origin="upper")
    ax.set_title(f"PAE — {pair}")
    ax.set_xlabel("residue (aligned)")
    ax.set_ylabel("residue (scored)")
    if chain_breaks:
        for b in chain_breaks:
            ax.axhline(b - 0.5, color="black", lw=0.8)
            ax.axvline(b - 0.5, color="black", lw=0.8)
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("PAE (Å)")
    plt.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def main():
    QCDIR.mkdir(exist_ok=True)
    manifest = pd.read_csv(MANIFEST, sep="\t")
    manifest_lookup = {r.pair: r for r in manifest.itertuples(index=False)}

    records = []
    contacts_rows = []
    # Iterate over canonical pair dirs (no trailing _<timestamp>); skip AF3's sibling
    # output dirs which match <pair>_2YYYYMMDD_*.
    pair_dirs = [p for p in sorted(RESULTS.iterdir())
                 if p.is_dir() and p.name in manifest_lookup]
    for pair_dir in pair_dirs:
        pair = pair_dir.name
        rec = {"pair": pair}
        meta = manifest_lookup.get(pair)
        if meta:
            rec.update({"p1": meta.p1, "p2": meta.p2,
                        "len1": meta.len1, "len2": meta.len2,
                        "total_tokens": meta.total_tokens})

        rec["msa_done"] = (pair_dir / f"{pair}_data.json").exists()
        af3_out = find_af3_output(pair_dir)
        rec["inf_done"] = af3_out is not None

        summary = load_summary(af3_out, pair) if af3_out else None
        if summary:
            for k in ("iptm", "ptm", "ranking_score",
                      "fraction_disordered", "has_clash", "num_recycles"):
                rec[k] = summary.get(k)
            chain_iptm = summary.get("chain_iptm") or []
            for i, v in enumerate(chain_iptm):
                rec[f"chain_iptm_{i}"] = v

            full = load_full_confidences(af3_out, pair)
            cif = find_model_cif(af3_out, pair)
            if cif:
                m = interface_metrics(cif)
                rec.update({f"iface_{k}": v for k, v in m.items()})
                contacts_rows.append({"pair": pair, **m})
            if full and "pae" in full:
                pae = np.asarray(full["pae"])
                token_chain = full.get("token_chain_ids") or []
                breaks = []
                if token_chain:
                    for i in range(1, len(token_chain)):
                        if token_chain[i] != token_chain[i-1]:
                            breaks.append(i)
                plot_pae(pae, pair, QCDIR / f"pae_{pair}.png", chain_breaks=breaks)
        records.append(rec)

    df = pd.DataFrame(records)
    # Stable column order
    front = [c for c in ["pair", "p1", "p2", "len1", "len2", "total_tokens",
                          "msa_done", "inf_done",
                          "iptm", "ptm", "ranking_score",
                          "fraction_disordered", "has_clash",
                          "iface_n_contacts", "iface_min_dist"] if c in df.columns]
    rest = [c for c in df.columns if c not in front]
    df = df[front + rest]
    if "ranking_score" in df.columns:
        df = df.sort_values("ranking_score", ascending=False, na_position="last")
    df.to_csv(QCDIR / "qc_summary.tsv", sep="\t", index=False)
    print(f"Wrote {QCDIR/'qc_summary.tsv'} ({len(df)} pairs)", file=sys.stderr)

    if contacts_rows:
        pd.DataFrame(contacts_rows).to_csv(QCDIR / "interface_contacts.tsv",
                                            sep="\t", index=False)

    # Heatmaps
    proteins = sorted(set(list(manifest.p1) + list(manifest.p2)))
    for metric, fname, vmin, vmax, cmap in [
        ("iptm", "iptm_heatmap.pdf", 0, 1, "viridis"),
        ("ranking_score", "ranking_score_heatmap.pdf", 0, 1, "magma"),
    ]:
        if metric not in df.columns:
            continue
        n = len(proteins)
        grid = np.full((n, n), np.nan)
        for _, r in df.iterrows():
            v = r.get(metric)
            if pd.isna(v) or pd.isna(r.get("p1")) or pd.isna(r.get("p2")):
                continue
            i = proteins.index(r["p1"]); j = proteins.index(r["p2"])
            grid[i, j] = v; grid[j, i] = v
        fig, ax = plt.subplots(figsize=(9, 8))
        im = ax.imshow(grid, cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_xticks(range(n)); ax.set_yticks(range(n))
        ax.set_xticklabels(proteins, rotation=45, ha="right")
        ax.set_yticklabels(proteins)
        for i in range(n):
            for j in range(n):
                if not np.isnan(grid[i, j]):
                    ax.text(j, i, f"{grid[i,j]:.2f}", ha="center", va="center",
                            color="white" if grid[i, j] < (vmin+vmax)/2 else "black",
                            fontsize=7)
        ax.set_title(f"{metric} across all pairs (homodimers on diagonal)")
        cb = fig.colorbar(im, ax=ax); cb.set_label(metric)
        plt.tight_layout()
        fig.savefig(QCDIR / fname)
        plt.close(fig)
        print(f"Wrote {QCDIR/fname}", file=sys.stderr)

    # ipTM vs ranking_score scatter, outliers labelled
    scatter_out = QCDIR / "iptm_vs_ranking_scatter.pdf"
    plot_metric_scatter(df, scatter_out)
    print(f"Wrote {scatter_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
