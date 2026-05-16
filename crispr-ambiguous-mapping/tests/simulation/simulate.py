"""Simulator for CRISPR-Correct ground-truth FASTQs.

Generates a small paired-end FASTQ pair plus a truth dataframe describing every
read's source guide, edits, recombination state, UMI, and PCR copy index. The
FASTQ header layout matches the chrX-style real data so the existing barcode /
UMI regexes parse it correctly:

    header = f"@sim:{mol_id}:{copy}_{U6_FIXED}_{revcomp(guide_barcode)} 1:N:0:{UMI_6bp}AG+{I5_FIXED}"

- barcode regex  r"_([^_ ]+)[\\s+]"       captures revcomp(guide_barcode)
- UMI     regex  r":([^+:]{6})(.{2})\\+"  captures the 6 bp UMI

See ../README.md for full context.
"""
from __future__ import annotations

import random
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Fixed scaffold constants (real HbF-style library)
# ---------------------------------------------------------------------------

U6_FIXED = "TTT"                       # 3 bp token in the header; not in reads
I5_FIXED = "TCGAACTG"                  # 8 bp; constant per run (single sample)
UMI_LINKER = "AG"                      # 2 bp anchor after the 6 bp UMI

MID_SCAFFOLD = (
    "GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAG"
    "GCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTT"
)  # 88 bp — between protospacer (R1) and surrogate (R2 revcomp)

READ_LENGTH = 150
QUAL = "I" * READ_LENGTH               # Phred 40 constant

# Edit window (0-indexed, inclusive)
PROTO_EDIT_WINDOW = (2, 10)
SURR_EDIT_WINDOW = (8, 16)   # = protospacer window shifted by the 6 bp surrogate upstream

DNA_BASES = np.array(list("ACGT"))
RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(RC_TABLE)[::-1]


def _random_seq(length: int, rng: np.random.Generator) -> str:
    return "".join(rng.choice(DNA_BASES, size=length))


# ---------------------------------------------------------------------------
# Library generation
# ---------------------------------------------------------------------------

def _min_pairwise_hamming(existing: List[str], candidate: str) -> int:
    if not existing:
        return len(candidate)
    arr_e = np.frombuffer("".join(existing).encode(), dtype=np.uint8).reshape(len(existing), -1)
    arr_c = np.frombuffer(candidate.encode(), dtype=np.uint8)
    return int((arr_e != arr_c).sum(axis=1).min())


def _generate_unique_protospacers(n: int, rng: np.random.Generator, min_hamming: int = 4, length: int = 20) -> List[str]:
    chosen: List[str] = []
    attempts = 0
    while len(chosen) < n:
        attempts += 1
        if attempts > 100 * n:
            raise RuntimeError(f"Failed to generate {n} unique protospacers at min_hamming={min_hamming}")
        cand = _random_seq(length, rng)
        if _min_pairwise_hamming(chosen, cand) >= min_hamming:
            chosen.append(cand)
    return chosen


def _generate_unique_barcodes(n: int, rng: np.random.Generator, length: int = 4) -> List[str]:
    """4 bp has only 256 variants — we insist each barcode is unique but allow
    low pairwise hamming since with 45 guides we can't guarantee much more."""
    chosen = set()
    while len(chosen) < n:
        chosen.add(_random_seq(length, rng))
    return list(chosen)


def _derive_editable_igrna(source_proto: str, rng: np.random.Generator) -> Tuple[str, int, str]:
    """Pick an A position in [2,10], substitute to a random non-A base."""
    a_positions = [i for i in range(PROTO_EDIT_WINDOW[0], PROTO_EDIT_WINDOW[1] + 1) if source_proto[i] == "A"]
    if not a_positions:
        return None, -1, ""  # caller retries with different source
    pos = int(rng.choice(a_positions))
    new_base = str(rng.choice(np.array(["C", "G", "T"])))
    new_proto = source_proto[:pos] + new_base + source_proto[pos + 1 :]
    return new_proto, pos, new_base


def _derive_noneditable_igrna(source_proto: str, rng: np.random.Generator) -> Tuple[str, int, str]:
    """Two sub-types, 50/50: (a) non-A in window [2,10] or (b) any base outside window."""
    if rng.random() < 0.5:
        # (a) non-A inside window
        cands = [i for i in range(PROTO_EDIT_WINDOW[0], PROTO_EDIT_WINDOW[1] + 1) if source_proto[i] != "A"]
    else:
        # (b) any position outside window
        cands = [i for i in range(20) if not (PROTO_EDIT_WINDOW[0] <= i <= PROTO_EDIT_WINDOW[1])]
    if not cands:
        return None, -1, ""
    pos = int(rng.choice(cands))
    old_base = source_proto[pos]
    other_bases = [b for b in "ACGT" if b != old_base]
    new_base = str(rng.choice(np.array(other_bases)))
    new_proto = source_proto[:pos] + new_base + source_proto[pos + 1 :]
    return new_proto, pos, new_base


def _build_surrogate(protospacer: str, rng: np.random.Generator) -> str:
    """Surrogate = [6 bp upstream random][protospacer][6 bp PAM random]."""
    return _random_seq(6, rng) + protospacer + _random_seq(6, rng)


def generate_library(
    n_pgrna: int = 30,
    n_ig_editable: int = 10,
    n_ig_noneditable: int = 5,
    seed: int = 42,
) -> pd.DataFrame:
    """Return a DataFrame with columns protospacer, surrogate, barcode,
    guide_type (pgRNA/igRNA_editable/igRNA_noneditable), source_pgrna_id.

    Guides are ordered [pgRNA..., igRNA_editable..., igRNA_noneditable...].
    Barcodes are unique across the whole library; pgRNA and its derived igRNAs
    get distinct barcodes (enforced by uniqueness).
    """
    rng = np.random.default_rng(seed)
    total = n_pgrna + n_ig_editable + n_ig_noneditable

    # pgRNA protospacers
    pgrna_protos = _generate_unique_protospacers(n_pgrna, rng)

    # Barcodes: draw total unique 4 bp strings
    barcodes = _generate_unique_barcodes(total, rng)

    rows: List[dict] = []
    for i, proto in enumerate(pgrna_protos):
        rows.append({
            "protospacer": proto,
            "surrogate":   _build_surrogate(proto, rng),
            "barcode":     barcodes[i],
            "guide_type":  "pgRNA",
            "source_pgrna_id": i,
        })

    # editable igRNAs: each derived from a random pgRNA. Retry if source has no A in window.
    bc_i = n_pgrna
    for _ in range(n_ig_editable):
        for _attempt in range(200):
            src = int(rng.integers(0, n_pgrna))
            src_proto = rows[src]["protospacer"]
            new_proto, pos, new_base = _derive_editable_igrna(src_proto, rng)
            if new_proto is None:
                continue
            # Ensure uniqueness of protospacer across library
            if any(r["protospacer"] == new_proto for r in rows):
                continue
            rows.append({
                "protospacer": new_proto,
                "surrogate":   _build_surrogate(new_proto, rng),
                "barcode":     barcodes[bc_i],
                "guide_type":  "igRNA_editable",
                "source_pgrna_id": src,
            })
            bc_i += 1
            break
        else:
            raise RuntimeError("Failed to synthesize editable igRNA")

    # noneditable igRNAs
    for _ in range(n_ig_noneditable):
        for _attempt in range(200):
            src = int(rng.integers(0, n_pgrna))
            src_proto = rows[src]["protospacer"]
            new_proto, pos, new_base = _derive_noneditable_igrna(src_proto, rng)
            if new_proto is None:
                continue
            if any(r["protospacer"] == new_proto for r in rows):
                continue
            rows.append({
                "protospacer": new_proto,
                "surrogate":   _build_surrogate(new_proto, rng),
                "barcode":     barcodes[bc_i],
                "guide_type":  "igRNA_noneditable",
                "source_pgrna_id": src,
            })
            bc_i += 1
            break
        else:
            raise RuntimeError("Failed to synthesize noneditable igRNA")

    df = pd.DataFrame(rows)
    assert df["protospacer"].is_unique, "protospacer uniqueness violated"
    assert df["barcode"].is_unique, "barcode uniqueness violated"
    return df


# ---------------------------------------------------------------------------
# Read simulation
# ---------------------------------------------------------------------------

@dataclass
class SimConfig:
    n_molecules: int = 5000
    edit_prob_per_a: float = 0.15
    recomb_prob: float = 0.05
    igrna_assignment_prob: float = 0.2
    pcr_copy_dist: Tuple[Tuple[int, float], ...] = ((1, 0.5), (2, 0.3), (3, 0.2))
    seed: int = 42
    read_length: int = READ_LENGTH


def _apply_abe_edits(seq: str, window: Tuple[int, int], edit_prob: float, rng: np.random.Generator) -> Tuple[str, List[int]]:
    """A->G with probability edit_prob at each A in [window_lo, window_hi] inclusive."""
    seq_arr = list(seq)
    edits: List[int] = []
    for i in range(window[0], window[1] + 1):
        if seq_arr[i] == "A" and rng.random() < edit_prob:
            seq_arr[i] = "G"
            edits.append(i)
    return "".join(seq_arr), edits


def _make_header(mol_id: int, copy: int, guide_barcode: str, umi_6bp: str) -> str:
    return (
        f"@sim:{mol_id}:{copy}_{U6_FIXED}_{revcomp(guide_barcode)} "
        f"1:N:0:{umi_6bp}{UMI_LINKER}+{I5_FIXED}"
    )


def simulate_reads(library: pd.DataFrame, cfg: SimConfig) -> Tuple[pd.DataFrame, List[str], List[str]]:
    """Generate FASTQ lines + truth dataframe.

    Returns (truth_df, r1_lines, r2_lines) where each *_lines list already
    contains the 4-line-per-record formatted FASTQ text.
    """
    rng = np.random.default_rng(cfg.seed)

    # Useful views
    pgrna_idx = library.index[library["guide_type"] == "pgRNA"].tolist()
    igrna_idx = library.index[library["guide_type"] != "pgRNA"].tolist()
    all_idx = library.index.tolist()

    pcr_copies, pcr_weights = zip(*cfg.pcr_copy_dist)
    pcr_copies = np.asarray(pcr_copies)
    pcr_weights = np.asarray(pcr_weights) / sum(pcr_weights)

    r1_lines: List[str] = []
    r2_lines: List[str] = []
    truth_rows: List[dict] = []
    read_id = 0

    for mol_id in range(cfg.n_molecules):
        # Source guide
        if rng.random() < cfg.igrna_assignment_prob and igrna_idx:
            source_id = int(rng.choice(igrna_idx))
        else:
            source_id = int(rng.choice(pgrna_idx))
        src = library.iloc[source_id]

        # Apply independent edits to protospacer and surrogate
        edited_proto, proto_edits = _apply_abe_edits(src["protospacer"], PROTO_EDIT_WINDOW, cfg.edit_prob_per_a, rng)
        edited_surr, surr_edits = _apply_abe_edits(src["surrogate"], SURR_EDIT_WINDOW, cfg.edit_prob_per_a, rng)

        # Recombination: take surrogate + barcode from a different guide
        is_recomb = bool(rng.random() < cfg.recomb_prob)
        recomb_id: Optional[int] = None
        if is_recomb:
            pool = [i for i in all_idx if i != source_id]
            recomb_id = int(rng.choice(pool))
            partner = library.iloc[recomb_id]
            # Re-do surrogate edits against the partner's surrogate
            edited_surr, surr_edits = _apply_abe_edits(partner["surrogate"], SURR_EDIT_WINDOW, cfg.edit_prob_per_a, rng)
            observed_barcode = partner["barcode"]
        else:
            observed_barcode = src["barcode"]

        # UMI + PCR copies
        umi_6bp = _random_seq(6, rng)
        n_copies = int(rng.choice(pcr_copies, p=pcr_weights))

        # R1 body: protospacer + mid scaffold
        r1_body = (edited_proto + MID_SCAFFOLD)[: cfg.read_length]
        r1_body = r1_body.ljust(cfg.read_length, "N")
        # R2 body: revcomp(surrogate) + revcomp(mid_scaffold)
        r2_body = (revcomp(edited_surr) + revcomp(MID_SCAFFOLD))[: cfg.read_length]
        r2_body = r2_body.ljust(cfg.read_length, "N")

        for copy in range(n_copies):
            header = _make_header(mol_id, copy, observed_barcode, umi_6bp)

            r1_lines.append(header)
            r1_lines.append(r1_body)
            r1_lines.append("+")
            r1_lines.append(QUAL)

            r2_lines.append(header)
            r2_lines.append(r2_body)
            r2_lines.append("+")
            r2_lines.append(QUAL)

            truth_rows.append({
                "read_id": read_id,
                "molecule_id": mol_id,
                "pcr_copy": copy,
                "source_guide_id": source_id,
                "source_guide_type": src["guide_type"],
                "observed_protospacer": edited_proto,
                "observed_surrogate": edited_surr,
                "observed_guide_barcode": observed_barcode,
                "edit_positions_protospacer": proto_edits,
                "edit_positions_surrogate": surr_edits,
                "is_recombined": is_recomb,
                "recomb_partner_id": recomb_id,
                "umi": umi_6bp,
            })
            read_id += 1

    truth_df = pd.DataFrame(truth_rows)
    return truth_df, r1_lines, r2_lines


def write_fastqs(r1_lines: List[str], r2_lines: List[str], r1_path: str, r2_path: str) -> None:
    with open(r1_path, "w") as f:
        f.write("\n".join(r1_lines) + "\n")
    with open(r2_path, "w") as f:
        f.write("\n".join(r2_lines) + "\n")
