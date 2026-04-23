#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path


FOURFOLD_STARTS = {"GC", "CG", "GG", "CT", "CC", "TC", "AC", "GT"}


def read_fasta(path):
    seqs = {}
    order = []
    name = None
    chunks = []

    with open(path, "r") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                    order.append(name)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.replace(" ", "").upper())

        if name is not None:
            seqs[name] = "".join(chunks).upper()
            order.append(name)

    if not seqs:
        raise ValueError("No sequences found in FASTA")

    lengths = {len(seq) for seq in seqs.values()}
    if len(lengths) != 1:
        raise ValueError("All sequences in the MSA must have the same alignment length")

    return seqs, order


def first_nongap_base_index(seq):
    for i, b in enumerate(seq):
        if b != "-":
            return i
    return None


def first_nongap_codon(seq):
    start = first_nongap_base_index(seq)
    if start is None:
        return None, None
    codon = seq[start:start+3]
    return start, codon


def sanity_check_start_atg(seqs, min_atg=None):
    aln_len = len(next(iter(seqs.values())))
    names = list(seqs.keys())
    nseq = len(names)

    if min_atg is None:
        min_atg = nseq

    best_start = None

    for i in range(0, aln_len - 2):
        codons = [seqs[name][i:i+3] for name in names]
        atg_count = sum(1 for c in codons if c == "ATG")

        if atg_count >= min_atg:
            best_start = i
            break

    if best_start is None:
        raise ValueError(
            f"No position found where ≥{min_atg} sequences have ATG"
        )

    print(f"ATG start found at position {best_start+1} with {atg_count}/{nseq} sequences")

    return best_start


def sanity_check_divisible_by_three(aln_len, coding_start):
    coding_len = aln_len - coding_start
    if coding_len % 3 != 0:
        raise ValueError(
            f"Aligned coding span from ATG start is not divisible by 3: "
            f"alignment_length={aln_len}, coding_start_1based={coding_start + 1}, "
            f"coding_span={coding_len}"
        )


def is_clean_codon(codon):
    return len(codon) == 3 and all(base in "ACGT" for base in codon)


def scan_shared_fourfold(seqs, order, coding_start):
    aln_len = len(next(iter(seqs.values())))
    kept_nt_positions_1based = []
    kept_codon_indices_1based = []

    codon_idx = 0
    for i in range(coding_start, aln_len, 3):
        codon_idx += 1
        codons = [seqs[name][i:i+3] for name in order]

        if not all(is_clean_codon(c) for c in codons):
            continue

        starts = {c[:2] for c in codons}
        if len(starts) != 1:
            continue

        start2 = next(iter(starts))
        if start2 not in FOURFOLD_STARTS:
            continue

        kept_nt_positions_1based.append(i + 3)
        kept_codon_indices_1based.append(codon_idx)

    return kept_nt_positions_1based, kept_codon_indices_1based


def extract_sites(seqs, order, kept_nt_positions_1based):
    idxs = [p - 1 for p in kept_nt_positions_1based]
    out = {}
    for name in order:
        out[name] = "".join(seqs[name][i] for i in idxs)
    return out


def trim_to_shared_start(seqs, order, coding_start):
    trimmed = {}
    for name in order:
        trimmed[name] = seqs[name][coding_start:]
    return trimmed


def write_fasta(path, seqs, order, width=80):
    with open(path, "w") as out:
        for name in order:
            out.write(f">{name}\n")
            seq = seqs[name]
            for i in range(0, len(seq), width):
                out.write(seq[i:i+width] + "\n")


def normalize_output_path(out):
    p = Path(out)
    if p.suffix in [".fa", ".fasta"]:
        return str(p.with_suffix(".4fold.fa"))
    else:
        return str(p) + ".4fold.fa"


def write_lines(path, values):
    with open(path, "w") as out:
        for v in values:
            out.write(f"{v}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Extract shared 4-fold degenerate sites from a codon-aware MSA FASTA"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input codon-aware MSA in FASTA format"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output FASTA containing extracted shared 4-fold sites"
    )
    parser.add_argument(
        "--positions",
        default=None,
        help="Optional output file for retained nucleotide positions (1-based)"
    )
    parser.add_argument(
        "--codons",
        default=None,
        help="Optional output file for retained codon indices relative to ATG start (1-based)"
    )
    parser.add_argument(
    "-m", "--min_atg",
    type=int,
    default=None,
    help="Minimum number of sequences required to share ATG at a position (default: all)"
    )

    args = parser.parse_args()

    try:
        seqs, order = read_fasta(args.input)
        aln_len = len(next(iter(seqs.values())))

        coding_start = sanity_check_start_atg(seqs, args.min_atg)
        sanity_check_divisible_by_three(aln_len, coding_start)

        trimmed = trim_to_shared_start(seqs, order, coding_start)
        write_fasta(args.output + ".trimmed.fa", trimmed, order)

        kept_nt_positions, kept_codon_indices = scan_shared_fourfold(seqs, order, coding_start)
        extracted = extract_sites(seqs, order, kept_nt_positions)

        out_fa = normalize_output_path(args.output)
        write_fasta(out_fa, extracted, order)

        pos_path = args.positions
        if pos_path is None:
            pos_path = out_fa + ".positions.txt"

        codon_path = args.codons
        if codon_path is None:
            codon_path = out_fa + ".codons.txt"

        write_lines(pos_path, kept_nt_positions)
        write_lines(codon_path, kept_codon_indices)

        print(f"Input sequences: {len(order)}")
        print(f"Alignment length (nt): {aln_len}")
        print(f"Aligned ATG start (1-based): {coding_start + 1}")
        print(f"Shared 4-fold sites retained: {len(kept_nt_positions)}")
        print(f"Output FASTA: {out_fa}")
        print(f"Output nt positions: {pos_path}")
        print(f"Output codon indices: {codon_path}")

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()