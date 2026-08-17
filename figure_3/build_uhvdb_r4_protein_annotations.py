#!/usr/bin/env python
"""Build UHVDB r4 per-protein annotations TSV matching the r5 schema.

Annotation columns are rebuilt from the packaged r4 release tables (same logic as
toolkit/bin/uhvdb_build_metadata.py). Genomic coordinates (start/end/strand/partial)
are not present in the packaged release tables, so they are recovered from the r5
protein_annotations file for overlapping protein_ids.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import polars as pl

from build_uhvdb_r4_metadata import OUT_PATH as R4_METADATA, RELEASE

FIG3 = Path("/mmfs1/gscratch/pedslabs_hoffman/carsonjm/CFPhageome/repos/UHVDB/uhvdb-manuscript/figure_3")
R5_PROTEIN_ANNOTATIONS = Path(
    "/mmfs1/gscratch/pedslabs_hoffman/carsonjm/CFPhageome/repos/UHVDB/toolkit2/databases/uhvdb/5.0/uhvdb_protein_annotations.tsv.gz"
)
OUT_PATH = FIG3 / "uhvdb_r4_protein_annotations.tsv.gz"

PROTEIN_COLUMNS = [
    "protein_id", "hash", "genomovar_rep", "bakta_acc", "foldseek_acc", "ips_id",
    "card_acc", "vfdb_acc", "pharokka_annot", "pharokka_category", "phold_category",
    "phold_annot", "empathi_annot", "start", "end", "strand", "partial",
]

ANNOTATION_COMPARE_COLS = [
    "hash", "genomovar_rep", "bakta_acc", "foldseek_acc", "ips_id", "card_acc",
    "vfdb_acc", "pharokka_annot", "pharokka_category", "phold_category",
    "phold_annot", "empathi_annot", "start", "end", "strand", "partial",
]


def get_genomovar_reps() -> set[str]:
    print("Loading r4 genomovar reps from metadata...")
    reps = (
        pl.read_csv(
            R4_METADATA,
            separator="\t",
            columns=["uhvdb_id", "genomovar_rep"],
            null_values=["NA", ""],
        )
        .filter(pl.col("uhvdb_id") == pl.col("genomovar_rep"))
        ["uhvdb_id"]
        .unique()
        .to_list()
    )
    print(f"  genomovar_reps={len(reps):,}")
    return set(reps)


def build_protein_annotations(genomovar_reps: set[str]) -> pl.DataFrame:
    print("Loading proteinhash + annotation tables...")
    prothash_df = (
        pl.read_csv(RELEASE / "uhvdb_proteinhash.tsv.gz", separator="\t")
        .unique("protein_id")
        .with_columns(pl.col("protein_id").str.replace(r"_[^_]*$", "").alias("genomovar_rep"))
        .filter(pl.col("genomovar_rep").is_in(genomovar_reps))
    )
    print(f"  proteinhash rows for reps={prothash_df.height:,}")
    hash_df = prothash_df.select("hash").unique()

    def filter_hashes(df: pl.DataFrame, hash_col: str = "hash") -> pl.DataFrame:
        return df.join(hash_df, left_on=hash_col, right_on="hash", how="semi")

    bakta_df = filter_hashes(
        pl.read_csv(
            RELEASE / "uhvdb_bakta.tsv.gz",
            separator="\t",
            columns=["Locus Tag", "Accession"],
            new_columns=["hash", "bakta_acc"],
            null_values=["-"],
            skip_rows=5,
            ignore_errors=True,
        )
    ).unique("hash")

    foldseek_df = filter_hashes(
        pl.read_csv(
            RELEASE / "uhvdb_foldseek.tsv.gz",
            separator="\t",
            has_header=False,
            columns=["column_1", "column_2"],
            new_columns=["hash", "foldseek_acc"],
        )
    ).unique("hash")

    ips_df = (
        filter_hashes(
            pl.read_csv(
                RELEASE / "uhvdb_interproscan.tsv.gz",
                separator="\t",
                ignore_errors=True,
                has_header=False,
                columns=["column_1", "column_5"],
                new_columns=["hash", "ips_id"],
            )
        )
        .group_by("hash")
        .agg(pl.col("ips_id").cast(pl.String).str.join(",").alias("ips_id"))
    )

    try:
        card_df = filter_hashes(
            pl.read_csv(
                RELEASE / "uhvdb_card.tsv.gz",
                separator="\t",
                ignore_errors=True,
                has_header=False,
                columns=["column_1", "column_2"],
                new_columns=["hash", "card_acc"],
            )
        ).unique("hash")
    except Exception:
        card_df = pl.DataFrame(schema={"hash": pl.String, "card_acc": pl.String})

    try:
        vfdb_df = filter_hashes(
            pl.read_csv(
                RELEASE / "uhvdb_vfdb.tsv.gz",
                separator="\t",
                ignore_errors=True,
                has_header=False,
                columns=["column_1", "column_2"],
                new_columns=["hash", "vfdb_acc"],
            )
        ).unique("hash")
    except Exception:
        vfdb_df = pl.DataFrame(schema={"hash": pl.String, "vfdb_acc": pl.String})

    pharokka_df = filter_hashes(
        pl.read_csv(
            RELEASE / "uhvdb_pharokka.tsv.gz",
            separator="\t",
            ignore_errors=True,
            columns=["ID", "annot", "category"],
            new_columns=["hash", "pharokka_annot", "pharokka_category"],
        )
        .filter(pl.col("pharokka_annot") != "hypothetical protein")
    ).unique("hash")

    phold_df = filter_hashes(
        pl.read_csv(
            RELEASE / "uhvdb_phold.tsv.gz",
            separator="\t",
            ignore_errors=True,
            columns=["cds_id", "function", "product"],
            new_columns=["hash", "phold_category", "phold_annot"],
        )
        .filter(pl.col("phold_category") != "unknown function")
    ).unique("hash")

    empathi_df = filter_hashes(
        pl.read_csv(
            RELEASE / "uhvdb_empathi.csv.gz",
            ignore_errors=True,
            columns=["", "Annotation"],
            new_columns=["hash", "empathi_annot"],
        )
    ).unique("hash")

    print("Joining annotation tables (bakta inner join, matching build script)...")
    combined = (
        prothash_df
        .join(bakta_df, on="hash", how="inner")
        .join(foldseek_df, on="hash", how="left")
        .join(ips_df, on="hash", how="left")
        .join(card_df, on="hash", how="left")
        .join(vfdb_df, on="hash", how="left")
        .join(pharokka_df, on="hash", how="left")
        .join(phold_df, on="hash", how="left")
        .join(empathi_df, on="hash", how="left")
    )
    print(f"  annotated proteins={combined.height:,}")
    return combined


def add_coords_from_r5(annotations: pl.DataFrame) -> pl.DataFrame:
    print("Recovering start/end/strand/partial from r5 for r4-range proteins...")
    # Filter by UHVDB numeric ID rather than a 33M-member is_in (much faster)
    max_r4_id = (
        annotations
        .select(pl.col("genomovar_rep").str.extract(r"UHVDB-(\d+)", 1).cast(pl.Int64).max())
        .item()
    )
    print(f"  max r4 genomovar_rep id={max_r4_id}")
    coords = (
        pl.scan_csv(
            R5_PROTEIN_ANNOTATIONS,
            separator="\t",
            null_values=["", "NA"],
        )
        .select(["protein_id", "genomovar_rep", "start", "end", "strand", "partial"])
        .with_columns(
            pl.col("genomovar_rep").str.extract(r"UHVDB-(\d+)", 1).cast(pl.Int64).alias("nid")
        )
        .filter(pl.col("nid") <= max_r4_id)
        .select(["protein_id", "start", "end", "strand", "partial"])
        .unique("protein_id")
        .collect()
    )
    print(f"  r5 coord rows in r4 id range={coords.height:,}")
    out = annotations.join(coords, on="protein_id", how="left")
    n_with = out.filter(pl.col("start").is_not_null()).height
    print(f"  coords joined for {n_with:,} / {out.height:,} proteins")
    return out.select(PROTEIN_COLUMNS)


def write_annotations(df: pl.DataFrame) -> None:
    print(f"Writing {OUT_PATH}...")
    with gzip.open(OUT_PATH, "wb") as f:
        df.write_csv(f, separator="\t")
    print(f"  rows={df.height:,}  size_mb={OUT_PATH.stat().st_size / 1e6:.1f}")


def compare_to_r5(r4: pl.DataFrame) -> None:
    print("Comparing r4 protein annotations to r5...")
    max_r4_id = (
        r4.select(pl.col("genomovar_rep").str.extract(r"UHVDB-(\d+)", 1).cast(pl.Int64).max()).item()
    )
    r5 = (
        pl.scan_csv(
            R5_PROTEIN_ANNOTATIONS,
            separator="\t",
            null_values=["", "NA"],
        )
        .with_columns(
            pl.col("genomovar_rep").str.extract(r"UHVDB-(\d+)", 1).cast(pl.Int64).alias("nid")
        )
        .filter(pl.col("nid") <= max_r4_id)
        .drop("nid")
        .collect()
        .unique("protein_id")
    )
    joined = r4.join(r5, on="protein_id", how="inner", suffix="_r5")
    print(f"  r4 proteins={r4.height:,}  overlapping in r5={joined.height:,}")

    print("\n  Annotation column match rates:")
    match_exprs = []
    for col in ANNOTATION_COMPARE_COLS:
        r5_col = f"{col}_r5"
        if r5_col not in joined.columns:
            continue
        match_expr = pl.col(col).cast(pl.String).eq_missing(pl.col(r5_col).cast(pl.String))
        rate = joined.select(match_expr.mean()).item()
        n_mismatch = joined.select((~match_expr).sum()).item()
        print(f"    {col}: {rate:.4%} match ({n_mismatch:,} mismatches)")
        match_exprs.append(match_expr)

    if match_exprs:
        fully = match_exprs[0]
        for expr in match_exprs[1:]:
            fully = fully & expr
        n_full = joined.select(fully.sum()).item()
        print(
            f"\n  Fully matching rows: {n_full:,} / {joined.height:,} "
            f"({n_full / joined.height:.4%})"
        )

    for col in ["bakta_acc", "empathi_annot", "pharokka_annot", "start", "ips_id"]:
        r5_col = f"{col}_r5"
        if r5_col not in joined.columns:
            continue
        mism = joined.filter(
            ~pl.col(col).cast(pl.String).eq_missing(pl.col(r5_col).cast(pl.String))
        ).select(["protein_id", col, r5_col]).head(5)
        if mism.height:
            print(f"\n  Sample mismatches for {col}:")
            print(mism)


def main() -> None:
    reps = get_genomovar_reps()
    annotations = build_protein_annotations(reps)
    with_coords = add_coords_from_r5(annotations)
    write_annotations(with_coords)
    compare_to_r5(with_coords)


if __name__ == "__main__":
    main()
