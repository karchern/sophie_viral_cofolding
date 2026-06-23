#!/usr/bin/env python3
"""Build an AF3 input JSON for co-folding all proteins in a single FASTA.

Identical sequences in the FASTA are merged into AF3 homomer shorthand
(one entry with `"id": [list, of, chain, ids]`). This means the MSA is
computed once per unique sequence, and AF3's paired-MSA step knows the
chains are identical, which improves inter-chain alignment quality.

Usage: _make_multifold_json.py <fasta_path> <out_json_path> <name>
"""
import json
import string
import sys
from collections import OrderedDict


def parse_fasta(path):
    out = []  # list of (header_name, seq) in input order
    name, buf = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(buf)))
                name = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line)
        if name is not None:
            out.append((name, "".join(buf)))
    return out


def main():
    fasta_path, out_path, name = sys.argv[1], sys.argv[2], sys.argv[3]
    entries = parse_fasta(fasta_path)
    if len(entries) == 0:
        sys.exit(f"no sequences parsed from {fasta_path}")
    if len(entries) > 26:
        sys.exit(f"too many chains ({len(entries)}); AF3 chain ids are A-Z")

    chain_ids = list(string.ascii_uppercase[:len(entries)])
    seq_to_ids = OrderedDict()
    seq_to_names = OrderedDict()
    for cid, (hdr, seq) in zip(chain_ids, entries):
        seq_to_ids.setdefault(seq, []).append(cid)
        seq_to_names.setdefault(seq, []).append(hdr)

    sequences = []
    summary_lines = []
    for seq, ids in seq_to_ids.items():
        proto = {"id": ids if len(ids) > 1 else ids[0], "sequence": seq}
        sequences.append({"protein": proto})
        summary_lines.append(
            f"  chain(s) {','.join(ids)}: {seq_to_names[seq][0]} "
            f"({'x' + str(len(ids)) if len(ids) > 1 else 'single'}, "
            f"len={len(seq)})"
        )

    obj = {
        "name": name,
        "modelSeeds": [1],
        "sequences": sequences,
        "dialect": "alphafold3",
        "version": 1,
    }
    with open(out_path, "w") as f:
        json.dump(obj, f, indent=2)
    print(f"Wrote {out_path}", file=sys.stderr)
    print(f"  {len(entries)} chain(s) across {len(sequences)} unique sequence(s):",
          file=sys.stderr)
    for line in summary_lines:
        print(line, file=sys.stderr)


if __name__ == "__main__":
    main()
