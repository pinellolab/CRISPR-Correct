# CLAUDE.md — CRISPR-Correct

Guidance for Claude Code sessions working in this repository.

## What this package is

**CRISPR-Correct** (PyPI: `crispr-ambiguous-mapping`, Pinello Lab) is a Python package for mapping observed CRISPR guide-RNA reads from FASTQs against a known guide library using **per-base Hamming distance**. It was designed to tolerate imperfect matches caused by base-editor self-editing or sequencing errors, and it supports:

- Protospacer + surrogate + guide barcode + guide UMI + sample (cell) barcode
- Regex / flank / fixed-position extraction per component
- IUPAC-ambiguous whitelist bases
- UMI-collapsed and UMI-noncollapsed counting
- Three ambiguity strategies (`ambiguous_ignored`, `ambiguous_accepted`, `ambiguous_spread`)
- AnnData output with log-fold-change layer

**Does not handle indels.** For indel-bearing reads, use [CRISPR BEAN](https://github.com/pinellolab/CRISPR-BEAN) or [CRISPR SURF](https://github.com/pinellolab/CRISPR-SURF/).

## Which branch is live

`master` is tagged `0.0.41zenodo` at **v0.0.199**, but development has moved on. The branch in actual use across this workstation's projects is:

> ### `origin/multi-sample-support` — v0.0.236, 9 commits ahead of master
>
> Adds single-cell / multi-sample support (`sample_barcode_*` params), multi-file FASTQ input (`fastq_r1_fns: List[str]`), and **renames** master's `barcode_*`/`umi_*` params to `guide_barcode_*`/`guide_umi_*`. Not backward-compatible.

Default to `multi-sample-support` when answering user questions or writing new code. See `../BRANCHES.md` for the full per-branch catalog.

## Python constraint

`>=3.8, <3.12`. Install PyTorch is **not** required (this package is numpy/scipy-based; Millipede is a separate package).

## Repo layout

```
CRISPR-Correct/
├── README.md                         Original README (documents the MASTER API — stale)
├── Dockerfile                        Python 3.10-slim, pins crispr-ambiguous-mapping==0.0.156
├── 20241011_CRISPR_Correct_CD19Demo.ipynb   Canonical demo notebook (CD19 base-editor screen)
├── examples/                         Terra/Firecloud workflow JSON
├── notebooks/                        Additional demo notebooks
├── LICENSE.txt
└── crispr-ambiguous-mapping/         ── The actual installable package lives here ──
    ├── pyproject.toml                Poetry; name: crispr-ambiguous-mapping
    ├── dist/                         Pre-built wheels (0.0.211 → 0.0.236 on multi-sample-support)
    ├── tests/__init__.py             Empty — no test suite configured
    └── crispr_ambiguous_mapping/
        ├── mapping/                  Top-level entry points
        ├── parsing/                  FASTQ / header / UMI-tools parsers
        ├── processing/               Encoding, Hamming inference, counting, editing profiles
        ├── quality_control/          Per-tier QC tallies
        ├── models/                   Dataclasses & type aliases
        ├── postprocessing/           AnnData + LFC
        ├── visualization/            Mutation histograms, LFC scatter, trinucleotide signatures
        └── utility/                  Pickle I/O, BAM→FASTQ, editing helpers
```

## Primary entry points

- **`crispr_ambiguous_mapping.mapping.get_whitelist_reporter_counts_from_fastq`**
  Canonical entry point. On multi-sample-support it takes `fastq_r1_fns: List[str]` (note plural) plus per-component parsing + threshold parameters. Returns a `WhitelistReporterCountsResult`.

- **`crispr_ambiguous_mapping.mapping.get_whitelist_reporter_counts_from_umitools_output`** *(deprecated)*
  Still works for UMI-tools-preprocessed FASTQs (igRNAModelling uses this). File is annotated `# Deprecated` on multi-sample-support — prefer the standard entry point for new code.

## Standard post-processing pipeline

```python
match_set   = cam.processing.get_matchset_alleleseries(...)
mutations   = cam.processing.get_mutation_profile(...)
linked      = cam.processing.tally_linked_mutation_count_per_sequence(...)
cam.visualization.plot_mutation_count_histogram(linked.protospacer_total_mutation_counter, ...)
cam.utility.calculate_average_editing_frequency(linked.protospacer_total_mutation_counter)
cam.utility.save_or_load_pickle(...)
```

Recipe detail and a migration checklist are in `../USAGE.md`.

## Rename gotcha on multi-sample-support

The `CountInput` dataclass on multi-sample-support renamed `contains_barcode` → `contains_guide_barcode` and `contains_umi` → `contains_guide_umi`, but the post-processing functions still take kwargs named `contains_barcode` / `contains_umi`. When chaining:

```python
# Correct on multi-sample-support:
get_matchset_alleleseries(..., contains_barcode=result.count_input.contains_guide_barcode, ...)
```

Silently fails on master-era notebooks ported forward without the remap.

## Where real usage lives

Four in-tree projects drive this package. Open `../EXAMPLES.md` for the pointer list, then open the specific driver file the user is asking about:

1. `2023_12_BB_BisiTilingScreen/` — protospacer-only, single-read (master API).
2. `2024_06_BB_HbF_Minitiling/.../20260330_TerraPipelineReplication/` — consumer of Terra-produced `result_.pickle`.
3. `2024_06_BB_HbF_Minitiling/.../20260330_igRNAModelling/` — full triplet + UMI sensor screen (deprecated UMI-tools entry point).
4. `2025_10_BB_scCRISPRCorrect_IGVF/` — single-cell with `sample_barcode_*` (multi-sample-support API).

## Tests / CI

`tests/__init__.py` is empty; no pytest fixtures, no GitHub Actions. Validation is done via the demo notebook and against real project outputs.

## Do / don't

- **Do** default to `origin/multi-sample-support` for user questions.
- **Do** direct the user to `../USAGE.md` §7 for the parameter cheat sheet.
- **Do** warn users that porting old notebooks requires the rename table in `../USAGE.md` §6.
- **Do** respect the Python `>=3.8,<3.12` constraint.
- **Don't** rely on `README.md` inside this repo for current API — it documents master (v0.0.199).
- **Don't** commit generated `.whl` / `.tar.gz` files; those in `dist/` are already tracked on multi-sample-support.
- **Don't** edit source without explicit user request — the package is published.

## Documentation map

Four supporting documents live one directory up in `CRISPR-Correct-Folder/`:

- `../ARCHITECTURE.md` — Module map, seven-stage pipeline, data models, Hamming/ambiguity algorithm.
- `../USAGE.md` — Install, three canonical call recipes, post-processing pipeline, rename gotcha, parameter cheat sheet, common pitfalls.
- `../BRANCHES.md` — Per-branch catalog, divergence from master, recommended install.
- `../EXAMPLES.md` — The four in-tree projects with driver-file pointers and parameter shapes.
