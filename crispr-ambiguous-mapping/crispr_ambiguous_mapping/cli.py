"""Flat-arg CLI for CRISPR-Correct.

§4.5: command-line entry point so Terra / Nextflow / Snakemake pipelines can
invoke mapping without writing Python. Argument names mirror
`api.ParsingConfig` fields 1:1 — every field becomes a `--flag value` option
with `_` → `-`. No config-file indirection.

Install it via `pip install -e .` (registers the `crispr-correct` console
script through the pyproject `[tool.poetry.scripts]` entry) and run:

```
crispr-correct map \\
    --r1 R1.fq.gz \\
    [--r2 R2.fq.gz] \\
    --library library.tsv \\
    --out result.pickle \\
    --protospacer-start-position 0 --protospacer-length 20 \\
    --is-protospacer-r1/--no-is-protospacer-r1 \\
    --protospacer-hamming-threshold-strict 7 \\
    [--retain-inference-results] [--cores 4]
```

`--r1` / `--r2` can be repeated for multi-file input. Boolean flags use
click's `--flag/--no-flag` convention so users can negate explicitly.
"""
from __future__ import annotations

import pickle
import sys
from dataclasses import fields as _dc_fields
from typing import Any, Dict

try:
    import click
except ImportError as _e:  # pragma: no cover
    raise ImportError(
        "The `crispr-correct` CLI requires `click`. Install via "
        "`pip install click` (or reinstall the package with its pinned deps)."
    ) from _e


# ---------------------------------------------------------------------------
# Dynamic option generation: one click.option per ParsingConfig field.
# ---------------------------------------------------------------------------

# Fields we don't expose on the CLI (they aren't ParsingConfig fields).
_SKIP_CLI_FIELDS = set()


def _resolve_type(tp):
    """Unwrap Optional[...] / Union[X, None] to the non-None inner type."""
    import typing
    origin = typing.get_origin(tp)
    if origin is typing.Union:
        args = [a for a in typing.get_args(tp) if a is not type(None)]
        return args[0] if args else str
    return tp


def _field_to_click_option(name: str, tp):
    """Turn one `ParsingConfig` type-resolved field into a `click.option(...)` decorator."""
    flag = name.replace("_", "-")
    opt_name = f"--{flag}"
    tp = _resolve_type(tp)
    if tp is bool:
        return click.option(opt_name + "/--no-" + flag, name, default=None, help=f"bool flag for `{name}`.")
    if tp is int:
        return click.option(opt_name, name, type=int, default=None, help=f"int kwarg `{name}`.")
    if tp is float:
        return click.option(opt_name, name, type=float, default=None, help=f"float kwarg `{name}`.")
    return click.option(opt_name, name, type=str, default=None, help=f"str kwarg `{name}`.")


def _apply_parsingconfig_options(cmd):
    """Decorator-stack: attach one click.option per ParsingConfig field to `cmd`."""
    import typing
    from .api import ParsingConfig
    hints = typing.get_type_hints(ParsingConfig)
    for f in _dc_fields(ParsingConfig):
        if f.name in _SKIP_CLI_FIELDS:
            continue
        cmd = _field_to_click_option(f.name, hints.get(f.name, str))(cmd)
    return cmd


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------


@click.group()
@click.version_option()
def main():
    """CRISPR-Correct — Hamming-distance-based CRISPR guide mapping."""


@main.command("map")
@click.option("--r1", "r1_fns", type=click.Path(exists=True, dir_okay=False), required=True, multiple=True, help="R1 FASTQ path(s); repeat for multi-file input.")
@click.option("--r2", "r2_fns", type=click.Path(exists=True, dir_okay=False), multiple=True, help="R2 FASTQ path(s); repeat. Omit for single-end.")
@click.option("--library", "library_path", type=click.Path(exists=True, dir_okay=False), required=True, help="Whitelist library TSV (columns: protospacer [surrogate] [barcode]).")
@click.option("--out", "out_path", type=click.Path(dir_okay=False), required=True, help="Output pickle path for the WhitelistReporterCountsResult.")
@_apply_parsingconfig_options
def map_cmd(r1_fns, r2_fns, library_path, out_path, **parsing_kwargs):
    """Map FASTQs against a whitelist library. See `--help` for all parsing kwargs."""
    import pandas as pd
    from .api import map_fastq, ParsingConfig

    # Drop unset (None) kwargs so ParsingConfig uses its dataclass defaults.
    cfg_kwargs: Dict[str, Any] = {k: v for k, v in parsing_kwargs.items() if v is not None}
    cfg = ParsingConfig(**cfg_kwargs)

    library_df = pd.read_csv(library_path, sep="\t")
    # Keep only the sequence columns we actually need — the package's
    # internal `strip_series` applies `.rstrip()` to every column and will
    # crash on non-string metadata columns (e.g. integer guide ids).
    _seq_cols = [c for c in ("protospacer", "guide_surrogate", "surrogate", "guide_barcode", "barcode") if c in library_df.columns]
    if _seq_cols:
        library_df = library_df[_seq_cols]
    click.echo(f"Loaded library with {len(library_df)} rows; using columns {list(library_df.columns)}.", err=True)

    result = map_fastq(
        whitelist_guide_reporter_df=library_df,
        fastq_r1_fns=list(r1_fns),
        fastq_r2_fns=list(r2_fns) if r2_fns else None,
        config=cfg,
    )

    with open(out_path, "wb") as fh:
        pickle.dump(result, fh)
    click.echo(f"Wrote result pickle to {out_path}.", err=True)


@main.command("count")
@click.option("--in", "in_path", type=click.Path(exists=True, dir_okay=False), required=True, help="Pickle produced by `crispr-correct map`.")
@click.option("--out", "out_path", type=click.Path(dir_okay=False), required=True, help="TSV destination (one row per tier x strategy).")
@click.option("--tier", "tier", type=str, default="protospacer_match_surrogate_match_barcode_match", show_default=True, help="MatchTier to extract.")
@click.option("--strategy", "strategy", type=str, default="ambiguous_accepted_umi_noncollapsed_counterseries", show_default=True, help="Ambiguity/UMI strategy attribute name.")
def count_cmd(in_path, out_path, tier, strategy):
    """Emit one tier's count Series as a 2-column TSV (index, count)."""
    import pandas as pd

    with open(in_path, "rb") as fh:
        result = pickle.load(fh)

    allw = result.all_match_set_whitelist_reporter_counter_series_results
    tier_wrap = getattr(allw, tier, None)
    if tier_wrap is None:
        raise click.ClickException(f"Tier `{tier}` not found on all_match_set_whitelist_reporter_counter_series_results.")
    series = getattr(tier_wrap, strategy, None)
    if series is None:
        raise click.ClickException(f"Strategy attribute `{strategy}` not found on tier `{tier}`.")

    if isinstance(series, pd.Series):
        df = series.reset_index()
    else:
        df = pd.DataFrame(list(series.items()), columns=["index", "count"])
    df.to_csv(out_path, sep="\t", index=False)
    click.echo(f"Wrote {len(df)} rows to {out_path}.", err=True)


@main.command("save")
@click.option("--in", "in_path", type=click.Path(exists=True, dir_okay=False), required=True, help="Source pickle produced by `crispr-correct map`.")
@click.option("--out-dir", "out_dir", type=click.Path(file_okay=False), required=True, help="Directory to write parquet + JSON artifacts.")
@click.option("--overwrite", is_flag=True, default=False, help="Overwrite non-empty --out-dir.")
def save_cmd(in_path, out_dir, overwrite):
    """§7.4: write a mapping result to a directory as parquet + JSON (cross-language durable)."""
    from .api import save as api_save
    with open(in_path, "rb") as fh:
        result = pickle.load(fh)
    api_save(result, out_dir, overwrite=overwrite)
    click.echo(f"Saved mapping-result directory to {out_dir}.", err=True)


@main.command("load")
@click.option("--in-dir", "in_dir", type=click.Path(exists=True, file_okay=False), required=True, help="Directory produced by `crispr-correct save`.")
@click.option("--out", "out_path", type=click.Path(dir_okay=False), required=True, help="Destination pickle path.")
def load_cmd(in_dir, out_path):
    """§7.4: reconstruct a mapping result pickle from a save-directory."""
    from .api import load as api_load
    result = api_load(in_dir)
    with open(out_path, "wb") as fh:
        pickle.dump(result, fh)
    click.echo(f"Loaded mapping-result directory {in_dir} -> {out_path}.", err=True)


@main.command("alleles")
@click.option("--in", "in_path", type=click.Path(exists=True, dir_okay=False), required=True, help="Pickle produced by `crispr-correct map --retain-inference-results`.")
@click.option("--out", "out_path", type=click.Path(dir_okay=False), required=True, help="Parquet destination for the first non-empty allele table.")
@click.option("--tier", "tier", type=str, default="protospacer_match_surrogate_match_barcode_match", show_default=True, help="MatchTier attribute to extract.")
@click.option("--ambiguity", "amb", type=click.Choice(["ignored", "accepted", "spread"]), default="accepted", show_default=True)
@click.option("--umi-strategy", "umi_strat", type=click.Choice(["", "collapsed", "noncollapsed"]), default="noncollapsed", show_default=True)
def alleles_cmd(in_path, out_path, tier, amb, umi_strat):
    """§4.5: run post-processing allele extraction on a retained mapping result."""
    import pandas as pd
    from .api import alleles as api_alleles
    with open(in_path, "rb") as fh:
        result = pickle.load(fh)
    ms = api_alleles(
        result,
        tier=tier,
        contains_guide_surrogate=result.count_input.contains_guide_surrogate,
        contains_guide_barcode=result.count_input.contains_guide_barcode,
        contains_guide_umi=result.count_input.contains_guide_umi,
    )
    strat_attr = "ambiguous_" + amb + ("_umi_" + umi_strat if umi_strat else "") + "_allele_df"
    df = getattr(ms, strat_attr, None)
    if df is None:
        raise click.ClickException(f"allele attribute `{strat_attr}` is None on the match-set result.")
    df.to_parquet(out_path)
    click.echo(f"Wrote {len(df)} allele rows to {out_path}.", err=True)


if __name__ == "__main__":
    main()
