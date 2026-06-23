#!/usr/bin/env python3
"""Generate AF3 input JSONs for all pairwise combinations (with replacement) of proteins.

Reads $VCO_ROOT/proteins.faa, writes $VCO_ROOT/data/pair_jsons/<pair>.json
and a tab-delimited $VCO_ROOT/data/pair_manifest.tsv with lengths + token counts
(used downstream by submit_all.sh for GPU routing).
"""
import json
import os
import re
import sys
from itertools import combinations_with_replacement
from pathlib import Path

ROOT = Path(os.environ.get(
    "VCO_ROOT",
    "/g/typas/Personal_Folders/Nic/sophie_viral_cofolding",
))
FASTA = ROOT / "proteins.faa"
OUTDIR = ROOT / "data" / "pair_jsons"


def sanitize(name: str) -> str:
    n = name.lower()
    n = re.sub(r"[^a-z0-9]+", "_", n).strip("_")
    return n


def parse_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    buf: list[str] = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line)
        if name is not None:
            seqs[name] = "".join(buf)
    return seqs


def main():
    print(f"VCO_ROOT={ROOT}", file=sys.stderr)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    seqs = parse_fasta(FASTA)
    names = sorted(seqs.keys())
    print(f"Loaded {len(names)} proteins (length): "
          + ", ".join(f"{n}({len(seqs[n])})" for n in names), file=sys.stderr)

    pairs = list(combinations_with_replacement(names, 2))
    print(f"Writing {len(pairs)} pair JSONs to {OUTDIR}", file=sys.stderr)

    manifest = []
    for p1, p2 in pairs:
        s1, s2 = sanitize(p1), sanitize(p2)
        pair_name = f"{s1}__{s2}"
        if p1 == p2:
            sequences = [{"protein": {"id": ["A", "B"], "sequence": seqs[p1]}}]
            total_tokens = 2 * len(seqs[p1])  # homodimer is 2 copies
        else:
            sequences = [
                {"protein": {"id": "A", "sequence": seqs[p1]}},
                {"protein": {"id": "B", "sequence": seqs[p2]}},
            ]
            total_tokens = len(seqs[p1]) + len(seqs[p2])
        data = {
            "name": pair_name,
            "modelSeeds": [1],
            "sequences": sequences,
            "dialect": "alphafold3",
            "version": 1,
        }
        out = OUTDIR / f"{pair_name}.json"
        out.write_text(json.dumps(data, indent=2))
        manifest.append((pair_name, p1, p2, len(seqs[p1]), len(seqs[p2]), total_tokens))

    manifest_path = ROOT / "data" / "pair_manifest.tsv"
    with open(manifest_path, "w") as f:
        f.write("pair\tp1\tp2\tlen1\tlen2\ttotal_tokens\n")
        for row in manifest:
            f.write("\t".join(map(str, row)) + "\n")
    print(f"Wrote manifest {manifest_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
