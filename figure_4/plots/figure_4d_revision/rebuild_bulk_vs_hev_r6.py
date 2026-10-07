#!/usr/bin/env python3
"""Rebuild figure_4d Bulk vs HEV selected violins from UHVDB r6 CoverM + gene_coverage."""
from __future__ import annotations

import gzip
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
from scipy import stats as scipy_stats

ROOT = Path("/mmfs1/gscratch/pedslabs_hoffman/carsonjm/CFPhageome/repos/UHVDB")
FIG_S15 = ROOT / "uhvdb-manuscript/figure_s15"
BULK_RUN = ROOT / "uhvdb-manuscript/figure_4/activity_profiling/bulk_figure2e_r6"
HEV_RUN = ROOT / "uhvdb-manuscript/figure_4/activity_profiling/other_hev_samples_r6"
METADATA = ROOT / "toolkit2/databases/uhvdb/v6/analyze/uhvdb_metadata.tsv.gz"
OUT = ROOT / "uhvdb-manuscript/figure_4/plots/figure_4d_revision"
OUT.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(FIG_S15))
from _he6_analysis_core import with_bio_group_flags  # noqa: E402

MIN_MEAN_COV = 5.0
GENE_BREADTH_THRESHOLD = 0.8
OUT_STEM = "figure_4d_selected_violins_bulk_vs_hev"

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
ZERO_COLS = [
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


def load_depth(path: Path, sample_id: str, arm: str) -> pl.DataFrame:
    open_fn = gzip.open if path.suffix == ".gz" or path.name.endswith(".gz") else open
    # polars reads gz natively
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
                pl.lit(arm).alias("arm"),
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


def list_samples(run: Path) -> list[str]:
    sheet = run / "refs" / "samples.tsv"
    lines = [ln for ln in sheet.read_text().splitlines() if ln.strip()]
    # bulk sheet has header; hev may not
    first = lines[0].split("\t")[0]
    start = 1 if first in {"sample", "sample_id"} else 0
    return [ln.split("\t")[0] for ln in lines[start:]]


def collect_depths(run: Path, arm: str) -> pl.DataFrame:
    parts = []
    n_ok = n_empty = n_miss = 0
    for sid in list_samples(run):
        depth = run / "results" / sid / f"{sid}.depth.tsv.gz"
        empty = run / "results" / sid / "no_contained_viruses"
        if depth.is_file():
            parts.append(load_depth(depth, sid, arm))
            n_ok += 1
        elif empty.exists():
            n_empty += 1
        else:
            n_miss += 1
    print(f"{arm}: depth={n_ok} empty={n_empty} missing={n_miss}")
    if not parts:
        raise SystemExit(f"No depths for {arm}")
    return pl.concat(parts, how="diagonal_relaxed")


def gene_coverage_counts(detections: pl.DataFrame, run_by_arm: dict[str, Path]) -> pl.DataFrame:
    need = detections.select(["sample_id", "contig_id", "arm"]).unique()
    pairs_by_sample: dict[str, tuple[str, set[str]]] = {}
    for g in need.partition_by("sample_id", as_dict=False):
        sid = g["sample_id"][0]
        arm = g["arm"][0]
        pairs_by_sample[sid] = (arm, set(g["contig_id"].to_list()))

    parts = []
    missing = 0
    for i, (sid, (arm, contig_ids)) in enumerate(sorted(pairs_by_sample.items()), 1):
        if i % 50 == 0 or i == 1:
            print(f"  gene_cov {i}/{len(pairs_by_sample)} {sid}", flush=True)
        path = run_by_arm[arm] / "results" / sid / f"{sid}.gene_coverage.tsv.gz"
        base = pl.DataFrame(
            {"sample_id": [sid] * len(contig_ids), "contig_id": list(contig_ids)}
        )
        if not path.is_file():
            missing += 1
            parts.append(base.with_columns([pl.lit(0).alias(c) for c in ZERO_COLS]))
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
            parts.append(base.with_columns([pl.lit(0).alias(c) for c in ZERO_COLS]))
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
            base.join(agg, on=["sample_id", "contig_id"], how="left").with_columns(
                [pl.col(c).fill_null(0) for c in ZERO_COLS]
            )
        )
    print(f"gene_cov missing files among detections: {missing}")
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


def main() -> None:
    meta = (
        pl.read_csv(
            METADATA,
            separator="\t",
            columns=[
                "uhvdb_id",
                "species_rep",
                "species_cluster_id",
                "ictv_class",
                "contig_length",
                "host_genes",
                "viral_genes",
                "checkv_quality",
                "aai_id",
                "aai_af",
            ],
            infer_schema_length=10_000,
        )
        .filter(pl.col("uhvdb_id") == pl.col("species_rep"))
        .filter(pl.col("ictv_class") == "Caudoviricetes")
        .with_columns(pl.col("species_cluster_id").cast(pl.Int64))
        .unique("uhvdb_id")
    )
    print("Caudoviricetes species-reps in metadata:", meta.height)

    bulk_depth = collect_depths(BULK_RUN, "bulk")
    hev_depth = collect_depths(HEV_RUN, "hev")
    coverm = pl.concat([bulk_depth, hev_depth], how="diagonal_relaxed")

    det = (
        coverm.filter(pl.col("contig_id").str.starts_with("UHVDB-"))
        .filter(pl.col("breadth").fill_null(0.0) > 0)
        .filter(pl.col("trimmed_mean").fill_null(0.0) >= MIN_MEAN_COV)
        .join(meta, left_on="contig_id", right_on="uhvdb_id", how="inner")
        .unique(["arm", "sample_id", "species_cluster_id"], keep="first")
        .with_columns(
            [
                pl.col("trimmed_mean").alias("mean"),
                pl.col("contig_length").cast(pl.Float64).alias("genome_length"),
                ((pl.col("aai_id") / 100.0) * pl.col("aai_af")).alias("aai_id_af"),
                (pl.col("checkv_quality") == "Complete")
                .cast(pl.Float64)
                .alias("complete_count"),
            ]
        )
    )
    print(
        det.group_by("arm").agg(
            [
                pl.len().alias("n_det"),
                pl.col("sample_id").n_unique().alias("n_samples"),
            ]
        )
    )

    gene = gene_coverage_counts(det, {"bulk": BULK_RUN, "hev": HEV_RUN})
    det = (
        det.join(gene, on=["sample_id", "contig_id"], how="left")
        .with_columns([pl.col(c).fill_null(0) for c in ZERO_COLS])
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
                    pl.col("host_genes").cast(pl.Float64)
                    * 10000.0
                    / pl.col("genome_length")
                ).alias("host_genes_per_10kb"),
                (
                    pl.col("viral_genes").cast(pl.Float64)
                    * 10000.0
                    / pl.col("genome_length")
                ).alias("viral_genes_per_10kb"),
            ]
        )
    )
    frac = (det["n_genes_breadth_ge80"] > 0).mean()
    print(f"frac rows with ≥1 gene breadth≥0.8: {frac:.3f}")

    paired_path = OUT / "bulk_vs_hev_ge5x_r6_genecov_ge80.tsv"
    keep = [
        "arm",
        "sample_id",
        "contig_id",
        "species_cluster_id",
        "trimmed_mean",
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
        "genome_length",
    ]
    det.select([c for c in keep if c in det.columns]).write_csv(paired_path, separator="\t")
    print("wrote", paired_path)

    left = det.filter(pl.col("arm") == "bulk")
    right = det.filter(pl.col("arm") == "hev")
    pal = {"left": "#4C72B0", "right": "#DD8452"}

    plt.rcParams.update(
        {
            "font.size": 12,
            "axes.labelsize": 12,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
        }
    )

    rows = []
    pvals = []
    for metric, label in SELECTED:
        a = np.asarray(left[metric].to_numpy(), dtype=float)
        b = np.asarray(right[metric].to_numpy(), dtype=float)
        a = a[np.isfinite(a)]
        b = b[np.isfinite(b)]
        p = (
            float(scipy_stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
            if a.size and b.size
            else float("nan")
        )
        pvals.append(p)
        rows.append(
            {
                "metric": metric,
                "label": label,
                "n_left": int(a.size),
                "n_right": int(b.size),
                "mannwhitney_p": p,
            }
        )
    qvals = bh_qvalues(pvals)
    for i, q in enumerate(qvals):
        rows[i]["bh_q"] = float(q)
        rows[i]["sig"] = sig_label(q)
    stats = pd.DataFrame(rows)
    stats_path = OUT / f"{OUT_STEM}_stats.tsv"
    stats.to_csv(stats_path, sep="\t", index=False)
    print(stats[["label", "n_left", "n_right", "bh_q", "sig"]].to_string(index=False))
    print("wrote", stats_path)

    n = len(SELECTED)
    ncol, nrow = 3, int(np.ceil(n / 3))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.4 * ncol, 3.6 * nrow), squeeze=False)
    qmap = dict(zip(stats["metric"], stats["bh_q"]))
    for ax, (metric, label) in zip(axes.ravel(), SELECTED):
        a = np.asarray(left[metric].to_numpy(), dtype=float)
        b = np.asarray(right[metric].to_numpy(), dtype=float)
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
        for body, key in zip(parts["bodies"], ["left", "right"]):
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
            sig_label(qmap[metric]),
            ha="center",
            va="bottom",
            fontsize=12,
        )
        ax.set_ylim(y0, y1 + 0.14 * span)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Bulk", "HEV"])
        ax.set_ylabel(label)
        ax.grid(axis="y", alpha=0.3)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    fig.tight_layout()
    png = OUT / f"{OUT_STEM}.png"
    pdf = OUT / f"{OUT_STEM}.pdf"
    fig.savefig(png, dpi=200, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print("wrote", png)
    print("wrote", pdf)

    n_b_samp = left["sample_id"].n_unique()
    n_h_samp = right["sample_id"].n_unique()
    note = OUT / f"{OUT_STEM}_NOTE.txt"
    note.write_text(
        "Bulk vs HEV selected violins rebuilt entirely from UHVDB r6:\n"
        f"- Bulk: figure2e bulk r6 CoverM trimmed_mean ≥{MIN_MEAN_COV:g}× "
        f"Caudoviricetes species-reps ({left.height} detections / {n_b_samp} samples).\n"
        f"- HEV: other_hev_samples_r6 CoverM ≥{MIN_MEAN_COV:g}× "
        f"Caudoviricetes ({right.height} detections / {n_h_samp} samples with depth).\n"
        "- Gene-category /10 kb + MCP+TerL+portal from r6 gene_coverage "
        f"(breadth ≥ {GENE_BREADTH_THRESHOLD}).\n"
        "- Host/viral /10 kb and AAI×AF from v6 CheckV metadata; breadth from r6 CoverM.\n"
        f"- Fraction of detection rows with ≥1 gene at breadth≥0.8: {frac:.3f}.\n"
    )
    print("wrote", note)


if __name__ == "__main__":
    main()
