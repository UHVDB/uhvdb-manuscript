#!/usr/bin/env python3
"""Rebuild figure_4d MMC PTH-FC≥1.5 selected violins entirely from UHVDB r6.

Uses:
  - r6 CoverM depths (virus ≥5× detections)
  - r6 sylph profiles (cascade PTH / fold-change induction)
  - r6 gene_coverage (breadth ≥ 0.8 gene-category metrics)
  - v6 UHVDB metadata (CheckV host/viral genes, AAI, ICTV, hosts)
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from scipy import stats as scipy_stats

BASE = Path("/mmfs1/gscratch/pedslabs_hoffman/carsonjm/CFPhageome/repos/UHVDB")
ANALYSIS = BASE / "uhvdb-manuscript-update/mmc_activity_analysis"
MMC_R6 = ANALYSIS / "mmc_r6/results"
SYLPH_DIR = MMC_R6 / "spring_sylph"
COVERM_DIR = MMC_R6 / "spring_coverm"
GENECOVER_DIR = MMC_R6 / "uhvdb_genecoverage"
SHEET = ANALYSIS / "paired_mmc_isolates_samplesheet.csv"
METADATA = BASE / "toolkit2/databases/uhvdb/v6/analyze/uhvdb_metadata.tsv.gz"
GTDB = BASE / "uhvdb-manuscript/figure_4/sylph_tax/gtdb_r226_metadata.tsv.gz"
FIG_S15 = BASE / "uhvdb-manuscript/figure_s15"
OUT = BASE / "uhvdb-manuscript/figure_4/plots/figure_4d_revision"
OUT.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(FIG_S15))
from _he6_analysis_core import with_bio_group_flags  # noqa: E402

MIN_MEAN_COV = 5.0
PTH_FC_THRESHOLD = float(__import__("os").environ.get("PTH_FC_THRESHOLD", "1.5"))
PTH_SENTINEL = 2.0
GENE_BREADTH_THRESHOLD = 0.8
MIN_GENOMOVAR_FRAC = 0.5
# Output stem override, e.g. figure_4d_selected_violins_no_mmc_vs_mmc
OUT_STEM = __import__("os").environ.get(
    "MMC_FIG_STEM", "figure_4d_selected_violins_mmc_pth_fc1.5"
)
PAIRED_STEM = __import__("os").environ.get(
    "MMC_PAIRED_STEM", "mmc_paired_uninduced_vs_induced_pth_fc1.5_nommc_metrics"
)


def pair_and_status(sample_id: str) -> tuple[str | None, str]:
    if sample_id.endswith("_no_mmc"):
        return "no_mmc", sample_id[: -len("_no_mmc")]
    if sample_id.endswith("_mmc"):
        return "mmc", sample_id[: -len("_mmc")]
    return None, sample_id


def load_depth(path: Path, sample_id: str, group: str) -> pl.DataFrame:
    raw = pl.read_csv(path, separator="\t")
    cols = raw.columns
    rename = {
        cols[0]: "contig_id",
        cols[1]: "trimmed_mean",
        cols[2]: "mean",
        cols[3]: "variance",
        cols[4]: "covered_bases",
        cols[5]: "length",
    }
    return (
        raw.rename(rename)
        .with_columns(
            [
                pl.col("trimmed_mean").cast(pl.Float64),
                pl.col("mean").cast(pl.Float64),
                pl.col("variance").cast(pl.Float64),
                pl.col("covered_bases").cast(pl.Float64),
                pl.col("length").cast(pl.Float64),
                pl.lit(sample_id).alias("sample_id"),
                pl.lit(group).alias("group"),
            ]
        )
        .with_columns(
            [
                (pl.col("covered_bases") / pl.col("length")).alias("breadth"),
                (1 - math.e ** (-0.833 * pl.col("mean"))).alias("expected_breadth"),
            ]
        )
        .with_columns(
            pl.when(pl.col("expected_breadth") > 1e-6)
            .then(pl.col("breadth") / pl.col("expected_breadth"))
            .otherwise(None)
            .alias("breadth_ratio")
        )
    )


def pick_host_at_rank(
    df: pl.DataFrame,
    phist_tax: str,
    phist_conn: str,
    phist_agree: str,
    crispr_tax: str,
    crispr_conn: str,
    crispr_agree: str,
    rank: str,
) -> pl.DataFrame:
    return (
        df.select(
            [
                "uhvdb_id",
                pl.lit(rank).alias("rank"),
                pl.col(phist_tax).str.replace(r"^(s|g|f)__", "").alias("phist_tax"),
                pl.col(phist_conn).fill_null(0).alias("phist_conn"),
                pl.col(phist_agree).fill_null(0.0).alias("phist_agree"),
                pl.col(crispr_tax).str.replace(r"^(s|g|f)__", "").alias("crispr_tax"),
                pl.col(crispr_conn).fill_null(0).alias("crispr_conn"),
                pl.col(crispr_agree).fill_null(0.0).alias("crispr_agree"),
            ]
        )
        .with_columns(
            pl.when(pl.col("phist_tax") == pl.col("crispr_tax"))
            .then(pl.col("phist_tax"))
            .when(
                (pl.col("phist_conn") * pl.col("phist_agree"))
                >= (pl.col("crispr_conn") * pl.col("crispr_agree"))
            )
            .then(pl.col("phist_tax"))
            .otherwise(pl.col("crispr_tax"))
            .alias("final_taxonomy")
        )
        .select(["uhvdb_id", "rank", "final_taxonomy"])
    )


def consensus_host(df: pl.DataFrame, col: str) -> pl.DataFrame:
    n_all = df.group_by("species_cluster_id").agg(pl.len().alias("n_genomovars"))
    return (
        df.filter(pl.col(col).is_not_null())
        .group_by(["species_cluster_id", col])
        .agg(pl.len().alias("n"))
        .join(n_all, on="species_cluster_id", how="left")
        .with_columns((pl.col("n") / pl.col("n_genomovars")).alias("frac"))
        .filter(pl.col("frac") >= MIN_GENOMOVAR_FRAC)
        .sort(["species_cluster_id", "n", col], descending=[False, True, False])
        .unique("species_cluster_id", maintain_order=True)
        .select(["species_cluster_id", col])
    )


def compute_cascade_vhr(
    sample_ids: list[str],
    id_map: pl.DataFrame,
    host: pl.DataFrame,
    gtdb: pl.DataFrame,
) -> pl.DataFrame:
    virus_rows, bac_rows = [], []
    for sid in sample_ids:
        path = SYLPH_DIR / f"{sid}.profile.tsv"
        if not path.is_file():
            continue
        df = pl.read_csv(path, separator="\t")
        virus_rows.append(
            df.filter(pl.col("Contig_name").str.starts_with("UHVDB-"))
            .with_columns(
                [
                    pl.lit(sid).alias("sample_id"),
                    pl.col("Contig_name").alias("uhvdb_id"),
                    pl.col("Taxonomic_abundance").cast(pl.Float64).alias("virus_tax_abund"),
                ]
            )
            .join(id_map, on="uhvdb_id", how="inner")
            .group_by(["sample_id", "species_cluster_id"])
            .agg(pl.col("virus_tax_abund").max())
        )
        bac_rows.append(
            df.filter(~pl.col("Contig_name").str.starts_with("UHVDB-"))
            .with_columns(
                [
                    pl.lit(sid).alias("sample_id"),
                    pl.col("Genome_file")
                    .cast(pl.Utf8)
                    .str.extract(r"(GC[AF]_\d+\.\d+)", 1)
                    .alias("accession"),
                    pl.col("Taxonomic_abundance").cast(pl.Float64).alias("host_tax_abund"),
                ]
            )
            .filter(pl.col("accession").is_not_null())
            .join(
                gtdb.select(["accession", "genus", "species", "family"]),
                on="accession",
                how="left",
            )
            .filter(pl.col("genus").is_not_null())
            .group_by(["sample_id", "species", "genus", "family"])
            .agg(pl.col("host_tax_abund").sum())
        )
    if not virus_rows:
        return pl.DataFrame(
            schema={
                "sample_id": pl.Utf8,
                "species_cluster_id": pl.Int64,
                "vhr_cascade": pl.Float64,
                "host_match_method": pl.Utf8,
            }
        )
    viruses = (
        pl.concat(virus_rows)
        .join(
            host.select(
                ["species_cluster_id", "final_species", "final_genus", "final_family"]
            ),
            on="species_cluster_id",
            how="left",
        )
        .with_columns(
            [
                pl.col("final_genus").alias("host_genus"),
                pl.col("final_family").alias("host_family"),
            ]
        )
    )
    bac_sylph = pl.concat(bac_rows)

    species_hits = (
        viruses.filter(pl.col("final_species").is_not_null() & (pl.col("virus_tax_abund") > 0))
        .join(
            bac_sylph.select(["sample_id", "species", "host_tax_abund"]),
            left_on=["sample_id", "final_species"],
            right_on=["sample_id", "species"],
            how="inner",
        )
        .with_columns(pl.lit("species_codetected").alias("host_match_method"))
    )
    matched_keys = species_hits.select(["sample_id", "species_cluster_id"]).unique()
    genus_singleton = (
        bac_sylph.group_by(["sample_id", "genus"])
        .agg(
            [
                pl.len().alias("n_species"),
                pl.col("species").first().alias("species"),
                pl.col("host_tax_abund").first().alias("host_tax_abund"),
            ]
        )
        .filter(pl.col("n_species") == 1)
    )
    genus_hits = (
        viruses.join(matched_keys, on=["sample_id", "species_cluster_id"], how="anti")
        .filter(pl.col("host_genus").is_not_null() & (pl.col("virus_tax_abund") > 0))
        .join(
            genus_singleton.select(["sample_id", "genus", "species", "host_tax_abund"]),
            left_on=["sample_id", "host_genus"],
            right_on=["sample_id", "genus"],
            how="inner",
        )
        .with_columns(pl.lit("genus_singleton").alias("host_match_method"))
    )
    matched_keys2 = pl.concat(
        [matched_keys, genus_hits.select(["sample_id", "species_cluster_id"]).unique()]
    ).unique()
    family_singleton = (
        bac_sylph.group_by(["sample_id", "family"])
        .agg(
            [
                pl.len().alias("n_species"),
                pl.col("species").first().alias("species"),
                pl.col("host_tax_abund").first().alias("host_tax_abund"),
            ]
        )
        .filter(pl.col("n_species") == 1)
    )
    family_hits = (
        viruses.join(matched_keys2, on=["sample_id", "species_cluster_id"], how="anti")
        .filter(pl.col("host_family").is_not_null() & (pl.col("virus_tax_abund") > 0))
        .join(
            family_singleton.select(["sample_id", "family", "species", "host_tax_abund"]),
            left_on=["sample_id", "host_family"],
            right_on=["sample_id", "family"],
            how="inner",
        )
        .with_columns(pl.lit("family_singleton").alias("host_match_method"))
    )
    vhr_cols = [
        "sample_id",
        "species_cluster_id",
        "virus_tax_abund",
        "host_tax_abund",
        "host_match_method",
    ]
    return (
        pl.concat(
            [
                species_hits.select(vhr_cols),
                genus_hits.select(vhr_cols),
                family_hits.select(vhr_cols),
            ]
        )
        .with_columns(
            (pl.col("virus_tax_abund") / pl.col("host_tax_abund")).alias("vhr_cascade")
        )
        .filter(pl.col("host_tax_abund") > 0)
        .sort(
            ["sample_id", "species_cluster_id", "virus_tax_abund"],
            descending=[False, False, True],
        )
        .unique(["sample_id", "species_cluster_id"], keep="first")
        .select(["sample_id", "species_cluster_id", "vhr_cascade", "host_match_method"])
    )


def gene_coverage_category_counts(detections: pl.DataFrame) -> pl.DataFrame:
    pairs_by_sample: dict[str, set[str]] = {}
    for g in detections.select(["sample_id", "contig_id"]).unique().partition_by(
        "sample_id", as_dict=False
    ):
        pairs_by_sample[g["sample_id"][0]] = set(g["contig_id"].to_list())

    zero_cols = [
        "n_genes_breadth_ge80",
        "num_capsid",
        "num_tail",
        "num_lysis",
        "n_integration",
        "n_amg_host_takeover",
        "n_dna_metabolism",
        "annot_mcp_hallmark",
        "annot_terl_hallmark",
        "annot_portal_hallmark",
        "n_conserved_hallmarks",
    ]
    parts = []
    missing = 0
    for sid, contig_ids in sorted(pairs_by_sample.items()):
        path = GENECOVER_DIR / f"{sid}.gene_coverage.tsv.gz"
        base_pairs = pl.DataFrame(
            {"sample_id": [sid] * len(contig_ids), "contig_id": list(contig_ids)}
        )
        if not path.is_file():
            missing += 1
            parts.append(base_pairs.with_columns([pl.lit(0).alias(c) for c in zero_cols]))
            continue
        gc = (
            pl.read_csv(
                path,
                separator="\t",
                columns=[
                    "genomovar_rep",
                    "pharokka_annot",
                    "pharokka_category",
                    "phold_category",
                    "empathi_annot",
                    "breadth",
                ],
            )
            .filter(
                (pl.col("breadth") >= GENE_BREADTH_THRESHOLD)
                & pl.col("genomovar_rep").is_in(list(contig_ids))
            )
            .with_columns(pl.lit(sid).alias("sample_id"))
        )
        if gc.height == 0:
            parts.append(base_pairs.with_columns([pl.lit(0).alias(c) for c in zero_cols]))
            continue
        flagged = with_bio_group_flags(gc)
        agg = (
            flagged.group_by(["sample_id", "genomovar_rep"])
            .agg(
                [
                    pl.len().alias("n_genes_breadth_ge80"),
                    pl.col("is_capsid_packaging").sum().cast(pl.Float64).alias("num_capsid"),
                    pl.col("is_tail").sum().cast(pl.Float64).alias("num_tail"),
                    pl.col("is_lysis").sum().cast(pl.Float64).alias("num_lysis"),
                    pl.col("is_integration").sum().cast(pl.Float64).alias("n_integration"),
                    pl.col("is_amg_host_takeover")
                    .sum()
                    .cast(pl.Float64)
                    .alias("n_amg_host_takeover"),
                    pl.col("is_dna_metabolism")
                    .sum()
                    .cast(pl.Float64)
                    .alias("n_dna_metabolism"),
                    (
                        ((pl.col("pharokka_annot") == "major head protein").sum() >= 1)
                        | ((pl.col("phold_category") == "major head protein").sum() >= 1)
                        | ((pl.col("empathi_annot") == "pvp|capsid|major_capsid").sum() >= 1)
                    )
                    .cast(pl.Int64)
                    .alias("annot_mcp_hallmark"),
                    (
                        ((pl.col("pharokka_annot") == "terminase large subunit").sum() >= 1)
                        | (
                            (pl.col("phold_category") == "terminase large subunit").sum()
                            >= 1
                        )
                        | (
                            (
                                pl.col("empathi_annot")
                                == "DNA-associated|terminase|packaging_assembly"
                            ).sum()
                            >= 1
                        )
                    )
                    .cast(pl.Int64)
                    .alias("annot_terl_hallmark"),
                    (
                        ((pl.col("pharokka_annot") == "portal protein").sum() >= 1)
                        | ((pl.col("phold_category") == "portal protein").sum() >= 1)
                        | ((pl.col("empathi_annot") == "pvp|portal").sum() >= 1)
                    )
                    .cast(pl.Int64)
                    .alias("annot_portal_hallmark"),
                ]
            )
            .rename({"genomovar_rep": "contig_id"})
            .with_columns(
                (
                    pl.col("annot_mcp_hallmark")
                    + pl.col("annot_terl_hallmark")
                    + pl.col("annot_portal_hallmark")
                ).alias("n_conserved_hallmarks")
            )
        )
        parts.append(
            base_pairs.join(agg, on=["sample_id", "contig_id"], how="left").with_columns(
                [pl.col(c).fill_null(0) for c in zero_cols]
            )
        )
    print(f"gene_cov: samples={len(pairs_by_sample)} missing_files={missing}")
    return pl.concat(parts, how="vertical_relaxed")


def bh_qvalues(pvals):
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    q_ranked = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = min(prev, ranked[i] * n / rank)
        q_ranked[i] = val
        prev = val
    out = np.empty(n)
    out[order] = q_ranked
    return np.clip(out, 0, 1)


def sig_label(q):
    if not np.isfinite(q):
        return "n/a"
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return "ns"


def cliffs_delta(x, y):
    """Cliff's δ for y − x (here: induced − uninduced)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if x.size == 0 or y.size == 0:
        return float("nan")
    diff = y[:, None] - x[None, :]
    return float((np.sum(diff > 0) - np.sum(diff < 0)) / (x.size * y.size))


def main() -> None:
    sheet = pl.read_csv(SHEET)
    samples = sheet["sample"].to_list()
    groups = sorted({pair_and_status(s)[1] for s in samples})

    complete = []
    for g in groups:
        s_no, s_mm = f"{g}_no_mmc", f"{g}_mmc"
        ok = all(
            (COVERM_DIR / f"{sid}.depth.tsv").is_file()
            and (SYLPH_DIR / f"{sid}.profile.tsv").is_file()
            and (GENECOVER_DIR / f"{sid}.gene_coverage.tsv.gz").is_file()
            for sid in (s_no, s_mm)
        )
        if ok:
            complete.append(g)
    print(
        f"sheet groups={len(groups)} complete r6 (coverm+sylph+genecov)={len(complete)} "
        f"missing={len(groups) - len(complete)}"
    )
    if not complete:
        raise SystemExit("No complete r6 pairs")

    # Metadata + hosts
    meta_cols = [
        "uhvdb_id",
        "genomovar_rep",
        "species_rep",
        "species_cluster_id",
        "ictv_class",
        "contig_length",
        "host_genes",
        "viral_genes",
        "checkv_quality",
        "virulent",
        "aai_id",
        "aai_af",
        "phist_gtdb_r226_species",
        "phist_species_connections",
        "phist_species_agreement",
        "phist_gtdb_r226_genus",
        "phist_genus_connections",
        "phist_genus_agreement",
        "phist_gtdb_r226_family",
        "phist_family_connections",
        "phist_family_agreement",
        "crispr_gtdb_r220_species",
        "crispr_species_connections",
        "crispr_species_agreement",
        "crispr_gtdb_r220_genus",
        "crispr_genus_connections",
        "crispr_genus_agreement",
        "crispr_gtdb_r220_family",
        "crispr_family_connections",
        "crispr_family_agreement",
    ]
    meta_full = pl.read_csv(METADATA, separator="\t", columns=meta_cols, infer_schema_length=10_000)
    id_map = (
        meta_full.select(["uhvdb_id", "species_cluster_id"])
        .unique("uhvdb_id")
        .with_columns(pl.col("species_cluster_id").cast(pl.Int64))
    )
    meta_sp = (
        meta_full.filter(pl.col("uhvdb_id") == pl.col("species_rep"))
        .select(
            [
                "uhvdb_id",
                "species_cluster_id",
                "ictv_class",
                "contig_length",
                "host_genes",
                "viral_genes",
                "checkv_quality",
                "virulent",
                "aai_id",
                "aai_af",
            ]
        )
        .with_columns(pl.col("species_cluster_id").cast(pl.Int64))
        .unique("uhvdb_id")
    )

    genomovar_meta = meta_full.filter(pl.col("uhvdb_id") == pl.col("genomovar_rep"))
    uhvdb_genomovar_host = (
        pl.concat(
            [
                pick_host_at_rank(
                    genomovar_meta,
                    "phist_gtdb_r226_species",
                    "phist_species_connections",
                    "phist_species_agreement",
                    "crispr_gtdb_r220_species",
                    "crispr_species_connections",
                    "crispr_species_agreement",
                    "species",
                ),
                pick_host_at_rank(
                    genomovar_meta,
                    "phist_gtdb_r226_genus",
                    "phist_genus_connections",
                    "phist_genus_agreement",
                    "crispr_gtdb_r220_genus",
                    "crispr_genus_connections",
                    "crispr_genus_agreement",
                    "genus",
                ),
                pick_host_at_rank(
                    genomovar_meta,
                    "phist_gtdb_r226_family",
                    "phist_family_connections",
                    "phist_family_agreement",
                    "crispr_gtdb_r220_family",
                    "crispr_family_connections",
                    "crispr_family_agreement",
                    "family",
                ),
            ]
        )
        .group_by("uhvdb_id")
        .agg(
            [
                pl.col("final_taxonomy")
                .filter(pl.col("rank") == "species")
                .first()
                .alias("final_species"),
                pl.col("final_taxonomy")
                .filter(pl.col("rank") == "genus")
                .first()
                .alias("final_genus"),
                pl.col("final_taxonomy")
                .filter(pl.col("rank") == "family")
                .first()
                .alias("final_family"),
            ]
        )
    )
    gv_in_species = (
        genomovar_meta.select(["uhvdb_id", "species_cluster_id", "species_rep"])
        .unique()
        .with_columns(pl.col("species_cluster_id").cast(pl.Int64))
        .join(uhvdb_genomovar_host, on="uhvdb_id", how="left")
        .rename({"uhvdb_id": "genomovar_rep"})
    )
    host = (
        gv_in_species.select(["species_cluster_id", "species_rep"])
        .unique()
        .join(consensus_host(gv_in_species, "final_species"), on="species_cluster_id", how="left")
        .join(consensus_host(gv_in_species, "final_genus"), on="species_cluster_id", how="left")
        .join(consensus_host(gv_in_species, "final_family"), on="species_cluster_id", how="left")
    )
    print(
        "host consensus:",
        host.filter(pl.col("final_species").is_not_null()).height,
        host.filter(pl.col("final_genus").is_not_null()).height,
        host.filter(pl.col("final_family").is_not_null()).height,
    )

    gtdb = (
        pl.read_csv(GTDB, separator="\t", has_header=False, new_columns=["accession", "lineage"])
        .with_columns(
            [
                pl.col("lineage").str.extract(r"g__([^;]+)", 1).alias("genus"),
                pl.col("lineage").str.extract(r"s__([^;]+)", 1).alias("species"),
                pl.col("lineage").str.extract(r"f__([^;]+)", 1).alias("family"),
            ]
        )
        .unique("accession")
    )

    # CoverM depths for complete pairs
    depth_lst = []
    for g in complete:
        for kind in ("no_mmc", "mmc"):
            sid = f"{g}_{kind}"
            depth_lst.append(load_depth(COVERM_DIR / f"{sid}.depth.tsv", sid, g))
    coverm_df = pl.concat(depth_lst)
    print("coverm rows", coverm_df.height, "samples", coverm_df["sample_id"].n_unique())

    # Species-rep Caudoviricetes detections
    det_base = (
        coverm_df.filter(pl.col("contig_id").str.starts_with("UHVDB-"))
        .filter(pl.col("breadth").fill_null(0.0) > 0)
        .join(meta_sp, left_on="contig_id", right_on="uhvdb_id", how="inner")
        .filter(pl.col("ictv_class") == "Caudoviricetes")
        .with_columns(pl.col("species_cluster_id").cast(pl.Int64))
        .unique(["group", "sample_id", "species_cluster_id"], keep="first")
    )
    no_mmc_all = det_base.filter(pl.col("sample_id").str.ends_with("_no_mmc"))
    mmc_all = det_base.filter(pl.col("sample_id").str.ends_with("_mmc"))
    print(
        f"Caudoviricetes detections: no_mmc={no_mmc_all.height} mmc={mmc_all.height}"
    )

    no_mmc_ge5 = no_mmc_all.filter(pl.col("trimmed_mean").fill_null(0.0) >= MIN_MEAN_COV)
    print(f"no_mmc ≥{MIN_MEAN_COV:g}×: {no_mmc_ge5.height} / {no_mmc_ge5['group'].n_unique()} groups")

    # Cascade PTH
    no_sids = [f"{g}_no_mmc" for g in complete]
    mm_sids = [f"{g}_mmc" for g in complete]
    vhr_no = compute_cascade_vhr(no_sids, id_map, host, gtdb)
    vhr_mm = compute_cascade_vhr(mm_sids, id_map, host, gtdb)
    print("cascade VHR rows: no_mmc", vhr_no.height, "mmc", vhr_mm.height)

    pth_no = vhr_no.select(
        [
            pl.col("sample_id").str.replace(r"_no_mmc$", "").alias("group"),
            "species_cluster_id",
            pl.col("vhr_cascade").alias("pth_no_mmc"),
        ]
    )
    pth_mm = vhr_mm.select(
        [
            pl.col("sample_id").str.replace(r"_mmc$", "").alias("group"),
            "species_cluster_id",
            pl.col("vhr_cascade").alias("pth_mmc"),
            pl.lit(True).alias("host_detected_mmc"),
        ]
    )
    mmc_virus_keys = (
        mmc_all.select(["group", "species_cluster_id"])
        .unique()
        .with_columns(pl.lit(True).alias("virus_detected_mmc"))
    )
    mmc_depth = (
        mmc_all.select(
            [
                "group",
                "species_cluster_id",
                pl.col("trimmed_mean").alias("trimmed_mean_mmc"),
                pl.col("breadth").alias("breadth_mmc"),
                pl.col("sample_id").alias("sample_id_mmc"),
            ]
        )
        .unique(["group", "species_cluster_id"], keep="first")
    )

    nommc = (
        no_mmc_ge5.join(pth_no, on=["group", "species_cluster_id"], how="left")
        .join(pth_mm, on=["group", "species_cluster_id"], how="left")
        .join(mmc_virus_keys, on=["group", "species_cluster_id"], how="left")
        .join(mmc_depth, on=["group", "species_cluster_id"], how="left")
        .with_columns(
            [
                pl.col("host_detected_mmc").fill_null(False),
                pl.col("virus_detected_mmc").fill_null(False),
            ]
        )
        .with_columns(
            # host undet sentinel: virus in MMC CoverM, no host-matched PTH
            (
                (~pl.col("host_detected_mmc")) & pl.col("virus_detected_mmc")
            ).alias("host_undetectable_mmc")
        )
        .with_columns(
            pl.when(pl.col("host_detected_mmc"))
            .then(pl.col("pth_mmc"))
            .when(pl.col("virus_detected_mmc"))
            .then(pl.lit(float(PTH_SENTINEL)))
            .otherwise(None)
            .alias("pth_mmc")
        )
        .with_columns((pl.col("pth_mmc") / pl.col("pth_no_mmc")).alias("pth_fold_change"))
        .with_columns(
            (
                pl.col("host_undetectable_mmc")
                | (
                    pl.col("pth_no_mmc").is_not_null()
                    & (pl.col("pth_no_mmc") > 0)
                    & (~pl.col("host_undetectable_mmc"))
                    & (pl.col("pth_fold_change") >= PTH_FC_THRESHOLD)
                )
            ).alias("induced")
        )
        .with_columns(
            pl.when(pl.col("induced"))
            .then(pl.lit("induced"))
            .otherwise(pl.lit("not_induced"))
            .alias("arm")
        )
        .with_columns(
            [
                pl.col("trimmed_mean").alias("mean"),
                ((pl.col("aai_id") / 100.0) * pl.col("aai_af")).alias("aai_id_af"),
                pl.col("contig_length").cast(pl.Float64).alias("genome_length"),
                (pl.col("checkv_quality") == "Complete").cast(pl.Float64).alias("complete_count"),
            ]
        )
    )

    print(
        nommc.select(
            [
                pl.len().alias("n"),
                (pl.col("arm") == "induced").sum().alias("n_induced"),
                (pl.col("arm") == "not_induced").sum().alias("n_uninduced"),
                pl.col("host_undetectable_mmc").sum().alias("n_host_undet"),
                (
                    (~pl.col("host_undetectable_mmc"))
                    & (pl.col("pth_fold_change") >= PTH_FC_THRESHOLD)
                )
                .sum()
                .alias("n_fc_ge"),
                pl.col("group").n_unique().alias("n_groups"),
            ]
        )
    )

    # Gene coverage metrics
    gene_counts = gene_coverage_category_counts(nommc)
    nommc = (
        nommc.join(gene_counts, on=["sample_id", "contig_id"], how="left")
        .with_columns(
            [
                pl.col(c).fill_null(0)
                for c in [
                    "n_genes_breadth_ge80",
                    "num_capsid",
                    "num_tail",
                    "num_lysis",
                    "n_integration",
                    "n_amg_host_takeover",
                    "n_dna_metabolism",
                    "n_conserved_hallmarks",
                ]
            ]
        )
        .with_columns(
            [
                (pl.col("num_capsid") * 10000.0 / pl.col("genome_length")).alias(
                    "capsid_per_10kb"
                ),
                (pl.col("num_tail") * 10000.0 / pl.col("genome_length")).alias(
                    "tail_per_10kb"
                ),
                (pl.col("num_lysis") * 10000.0 / pl.col("genome_length")).alias(
                    "lysis_per_10kb"
                ),
                (pl.col("n_integration") * 10000.0 / pl.col("genome_length")).alias(
                    "integration_per_10kb"
                ),
                (pl.col("n_amg_host_takeover") * 10000.0 / pl.col("genome_length")).alias(
                    "amg_host_takeover_per_10kb"
                ),
                (pl.col("n_dna_metabolism") * 10000.0 / pl.col("genome_length")).alias(
                    "dna_metabolism_per_10kb"
                ),
                (
                    pl.col("host_genes").cast(pl.Float64) * 10000.0 / pl.col("genome_length")
                ).alias("host_genes_per_10kb"),
                (
                    pl.col("viral_genes").cast(pl.Float64) * 10000.0 / pl.col("genome_length")
                ).alias("viral_genes_per_10kb"),
            ]
        )
    )

    frac_ge80 = (nommc["n_genes_breadth_ge80"] > 0).mean()
    print(f"frac rows with ≥1 gene breadth≥0.8: {frac_ge80:.3f}")

    keep_cols = [
        "group",
        "contig_id",
        "sample_id",
        "sample_id_mmc",
        "species_cluster_id",
        "breadth",
        "aai_id_af",
        "host_genes_per_10kb",
        "viral_genes_per_10kb",
        "lysis_per_10kb",
        "capsid_per_10kb",
        "tail_per_10kb",
        "amg_host_takeover_per_10kb",
        "dna_metabolism_per_10kb",
        "integration_per_10kb",
        "n_conserved_hallmarks",
        "n_genes_breadth_ge80",
        "trimmed_mean",
        "trimmed_mean_mmc",
        "breadth_mmc",
        "pth_mmc",
        "pth_no_mmc",
        "pth_fold_change",
        "host_undetectable_mmc",
        "host_detected_mmc",
        "virus_detected_mmc",
        "arm",
    ]
    # ensure sample_id_mmc filled
    nommc = nommc.with_columns(
        pl.col("sample_id_mmc").fill_null(pl.col("group") + "_mmc")
    )
    paired_path = OUT / f"{PAIRED_STEM}.tsv"
    nommc.select([c for c in keep_cols if c in nommc.columns]).write_csv(
        paired_path, separator="\t"
    )
    print("wrote", paired_path)

    # Plot
    SELECTED = [
        ("breadth", "Breadth"),
        ("aai_id_af", "CheckV AAI × AF to nearest DTR"),
        ("host_genes_per_10kb", "CheckV Host genes / 10 Kb"),
        ("viral_genes_per_10kb", "CheckV viral genes / 10 Kb"),
        ("lysis_per_10kb", "Lysis genes / 10 Kb"),
        ("capsid_per_10kb", "Capsid genes / 10 Kb"),
        ("tail_per_10kb", "Tail genes / 10 Kb"),
        ("amg_host_takeover_per_10kb", "AMG genes / 10 Kb"),
        ("dna_metabolism_per_10kb", "DNA metabolism genes / 10 Kb"),
        ("integration_per_10kb", "Integration genes / 10 Kb"),
        ("n_conserved_hallmarks", "MCP + TerL + portal"),
    ]
    LEFT, RIGHT = "not_induced", "induced"
    LEFT_LABEL, RIGHT_LABEL = "Uninduced", "Induced"
    pal = {LEFT: "#4C72B0", RIGHT: "#DD8452"}
    pdf = nommc.to_pandas()
    a_df = pdf[pdf["arm"] == LEFT]
    b_df = pdf[pdf["arm"] == RIGHT]

    plt.rcParams.update(
        {
            "font.size": 13,
            "axes.labelsize": 13,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
        }
    )
    SIG_FS = 15

    stats_rows = []
    pvals = []
    for col, _label in SELECTED:
        a = a_df[col].to_numpy(dtype=float)
        b = b_df[col].to_numpy(dtype=float)
        a = a[np.isfinite(a)]
        b = b[np.isfinite(b)]
        if a.size >= 1 and b.size >= 1:
            p = float(scipy_stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
        else:
            p = float("nan")
        pvals.append(p)
        stats_rows.append(
            {
                "metric": col,
                "n_uninduced": int(a.size),
                "n_induced": int(b.size),
                "median_uninduced": float(np.median(a)) if a.size else float("nan"),
                "median_induced": float(np.median(b)) if b.size else float("nan"),
                "mean_uninduced": float(np.mean(a)) if a.size else float("nan"),
                "mean_induced": float(np.mean(b)) if b.size else float("nan"),
                "p_raw": p,
                "cliffs_delta": cliffs_delta(a, b),  # induced − uninduced
            }
        )
    qvals = bh_qvalues(pvals)
    for i, q in enumerate(qvals):
        stats_rows[i]["q_bh"] = float(q)
        stats_rows[i]["sig"] = sig_label(q)

    stats_out = pl.DataFrame(
        [
            {
                "metric": r["metric"],
                "label": dict(SELECTED)[r["metric"]],
                "n_uninduced": r["n_uninduced"],
                "n_induced": r["n_induced"],
                "mannwhitney_p": r["p_raw"],
                "cliffs_delta_induced_minus_uninduced": r["cliffs_delta"],
                "bh_q": r["q_bh"],
                "sig": r["sig"],
            }
            for r in stats_rows
        ]
    )
    stats_path = OUT / f"{OUT_STEM}_stats.tsv"
    stats_out.write_csv(stats_path, separator="\t")
    print("wrote", stats_path)

    n = len(SELECTED)
    ncol, nrow = 3, int(np.ceil(n / 3))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 3.8 * nrow), squeeze=False)
    q_by_metric = {r["metric"]: r["q_bh"] for r in stats_rows}

    for ax, (metric, label) in zip(axes.ravel(), SELECTED):
        a = a_df[metric].to_numpy(dtype=float)
        b = b_df[metric].to_numpy(dtype=float)
        a = a[np.isfinite(a)]
        b = b[np.isfinite(b)]
        data = [a, b]
        parts = ax.violinplot(
            data,
            positions=[0, 1],
            showmeans=False,
            showmedians=False,
            showextrema=False,
            widths=0.9,
        )
        for body, key in zip(parts["bodies"], [LEFT, RIGHT]):
            body.set_facecolor(pal[key])
            body.set_edgecolor("black")
            body.set_alpha(0.55)
            body.set_linewidth(0.8)
        ax.boxplot(
            data,
            positions=[0, 1],
            widths=0.18,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "0.1", "linewidth": 1.3},
            whiskerprops={"color": "0.1", "linewidth": 1.0},
            capprops={"color": "0.1", "linewidth": 1.0},
            boxprops={
                "facecolor": "white",
                "edgecolor": "0.1",
                "linewidth": 1.0,
                "alpha": 0.9,
            },
        )
        ax.scatter(
            [0, 1],
            [float(np.mean(x)) for x in data],
            marker="D",
            s=28,
            color="0.1",
            zorder=5,
        )
        y0, y1 = ax.get_ylim()
        span = y1 - y0 if y1 > y0 else 1.0
        y_bar = y1 + 0.04 * span
        ax.plot(
            [0, 0, 1, 1],
            [y_bar - 0.015 * span, y_bar, y_bar, y_bar - 0.015 * span],
            color="0.1",
            lw=1.1,
        )
        ax.text(
            0.5,
            y_bar + 0.01 * span,
            sig_label(q_by_metric[metric]),
            ha="center",
            va="bottom",
            fontsize=SIG_FS,
        )
        ax.set_ylim(y0, y1 + 0.16 * span)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([LEFT_LABEL, RIGHT_LABEL])
        ax.set_ylabel(label)
        ax.tick_params(axis="both", labelsize=11)
        ax.yaxis.label.set_size(13)
        ax.grid(axis="y", alpha=0.3)

    for ax in axes.ravel()[n:]:
        ax.axis("off")

    fig.tight_layout()
    png = OUT / f"{OUT_STEM}.png"
    pdfp = OUT / f"{OUT_STEM}.pdf"
    fig.savefig(png, dpi=200, bbox_inches="tight")
    fig.savefig(pdfp, bbox_inches="tight")
    plt.close(fig)
    print("wrote", png, "figsize=", fig.get_size_inches())
    print("wrote", pdfp)

    n_u = int((pdf["arm"] == LEFT).sum())
    n_i = int((pdf["arm"] == RIGHT).sum())
    n_g = int(pdf["group"].nunique())
    n_undet = int(pdf["host_undetectable_mmc"].fillna(False).sum())
    n_fc = int(
        (
            (~pdf["host_undetectable_mmc"].fillna(False))
            & (pdf["pth_fold_change"] >= PTH_FC_THRESHOLD)
        ).sum()
    )
    note = OUT / f"{OUT_STEM}_NOTE.txt"
    note.write_text(
        f"Paired MMC selected violins with PTH fold-change induction (≥{PTH_FC_THRESHOLD:g}×), "
        "rebuilt entirely from UHVDB r6:\n"
        "- Unit: no_mmc ≥5× Caudoviricetes group×virus from r6 CoverM; "
        "metrics from the no_mmc sample.\n"
        "- Induction labels from r6 sylph cascade PTH: induced if\n"
        f"    (1) host detected in MMC and pth_mmc / pth_no_mmc ≥ {PTH_FC_THRESHOLD:g}, or\n"
        "    (2) host undetectable in MMC while virus is present in MMC CoverM "
        f"(sentinel pth_mmc={PTH_SENTINEL:g}).\n"
        "- Uninduced: remaining no_mmc ≥5× detections. MMC-only excluded.\n"
        "- Gene-category /10 kb + MCP+TerL+portal from no_mmc r6 gene_coverage "
        f"(breadth ≥ {GENE_BREADTH_THRESHOLD}).\n"
        "- Host/viral /10 kb and AAI×AF from v6 CheckV metadata; breadth from r6 CoverM.\n"
        f"- n: uninduced={n_u}, induced={n_i} "
        f"(FC≥{PTH_FC_THRESHOLD:g} host-detected={n_fc}; host-undet={n_undet}); groups={n_g}.\n"
        f"- Complete r6 pairs: {len(complete)}/94 "
        f"(missing: {sorted(set(groups) - set(complete))}).\n"
        "- Style/width matched to figure_4d_selected_violins panel.\n"
    )
    print("wrote", note)


if __name__ == "__main__":
    main()
