#!/usr/bin/env python3
"""Regenerate figure 2b/c genus and family accumulation curves for UHVDB r6."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import seaborn as sns
from matplotlib.ticker import MaxNLocator

BASE = Path(__file__).resolve().parent
META = (
    Path(__file__).resolve().parents[1]
    / "figure_1"
    / "uhvdb_v6_metadata_sra.tsv"
)
OUT = BASE / "plots"
OUT.mkdir(exist_ok=True)

SITE_ORDER = ["Skin", "Gut", "Urogenital", "Airways", "Blood"]
SITE_COLORS = {
    "Skin": "#1f77b4",
    "Gut": "#ff7f0e",
    "Urogenital": "#2ca02c",
    "Airways": "#d62728",
    "Blood": "#9467bd",
}


def load_unique_genomes() -> pl.DataFrame:
    return (
        pl.read_csv(
            META,
            separator="\t",
            columns=[
                "seq_name",
                "seqhash_rep",
                "body_site",
                "species_cluster_id",
                "genus_cluster_id",
                "family_cluster_id",
                "completeness_method",
            ],
            infer_schema_length=5000,
        )
        .filter(pl.col("seq_name") == pl.col("seqhash_rep"))
        .filter(pl.col("species_cluster_id").is_not_null())
        .with_columns(pl.col("completeness_method").fill_null(""))
    )


def load_species_table(genomes: pl.DataFrame) -> pl.DataFrame:
    """One row per species with genus/family + whether any member has TR."""
    return (
        genomes.group_by("species_cluster_id")
        .agg(
            [
                pl.col("genus_cluster_id").drop_nulls().first().alias("genus_cluster_id"),
                pl.col("family_cluster_id").drop_nulls().first().alias("family_cluster_id"),
                pl.col("completeness_method").str.contains("TR").any().alias("is_tr"),
            ]
        )
        .filter(pl.col("genus_cluster_id").is_not_null())
        .filter(pl.col("family_cluster_id").is_not_null())
    )


def accumulate_once(
    labels: np.ndarray, is_tr: np.ndarray, rng: np.random.Generator
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    import pandas as pd

    order = rng.permutation(labels.shape[0])
    s = pd.Series(labels[order])
    t = is_tr[order]
    total = (~s.duplicated()).cumsum().to_numpy(dtype=np.int32)
    non_sing = (s.groupby(s, sort=False).cumcount() == 1).cumsum().to_numpy(dtype=np.int32)
    s_tr = s.where(t)
    tr = (s_tr.notna() & ~s_tr.duplicated()).cumsum().to_numpy(dtype=np.int32)
    return total, non_sing, tr


def calculate_accumulation(
    df: pl.DataFrame, label_col: str, n_permutations: int = 50
) -> pl.DataFrame:
    labels = df[label_col].to_numpy()
    is_tr = df["is_tr"].to_numpy()
    n = labels.shape[0]
    sum_total = np.zeros(n, dtype=np.float64)
    sum_non = np.zeros(n, dtype=np.float64)
    sum_tr = np.zeros(n, dtype=np.float64)
    for i in range(n_permutations):
        total, non_sing, tr = accumulate_once(labels, is_tr, np.random.default_rng(i))
        sum_total += total
        sum_non += non_sing
        sum_tr += tr
        if (i + 1) % 10 == 0:
            print(f"  permutation {i + 1}/{n_permutations}")
    return pl.DataFrame(
        {
            "units_sampled": np.arange(1, n + 1),
            "mean_total": sum_total / n_permutations,
            "mean_non_singletons": sum_non / n_permutations,
            "mean_tr": sum_tr / n_permutations,
        }
    )


def plot_overall(
    curve: pl.DataFrame,
    *,
    ylabel: str,
    outfile: str,
    legend_loc: str = "upper left",
) -> None:
    sns.reset_orig()
    plt.rcParams.update({"font.size": 14})
    step = max(1, curve.height // 5000)
    plot_df = curve.gather_every(step)
    x = plot_df["units_sampled"].to_numpy()
    y_total = plot_df["mean_total"].to_numpy()
    y_non = plot_df["mean_non_singletons"].to_numpy()

    plt.figure(figsize=(7, 6))
    plt.grid(True, alpha=0.3)
    sns.lineplot(x=x, y=y_total, color="#2E86AB", linewidth=2.5, label="Total")
    sns.lineplot(x=x, y=y_non, color="#D62828", linewidth=2.5, label="Non-singleton")
    plt.xlabel("Number of species sampled", fontdict={"fontweight": "bold"})
    plt.ylabel(ylabel, fontdict={"fontweight": "bold"})
    plt.xlim(0, max(x))
    plt.ylim(0, max(y_total) * 1.05)
    plt.legend(title="Cluster type", loc=legend_loc)
    ax = plt.gca()
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
    plt.tight_layout()
    out = OUT / outfile
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"wrote {out}")


def body_site_cluster_accumulation(
    genomes: pl.DataFrame,
    cluster_col: str,
    *,
    n_reps: int = 50,
    n_points: int = 100,
) -> pl.DataFrame:
    site_df = (
        genomes.select(["seqhash_rep", "body_site", cluster_col])
        .filter(
            pl.col("body_site").is_not_null()
            & (pl.col("body_site") != "Other")
            & pl.col(cluster_col).is_not_null()
        )
        .unique()
    )
    results = []
    for body_site in site_df["body_site"].unique().to_list():
        d = site_df.filter(pl.col("body_site") == body_site)
        n_units = d.height
        if n_units == 0:
            continue
        labels = d[cluster_col].to_numpy()
        subset_sizes = np.unique(np.linspace(1, n_units, min(n_points, n_units), dtype=int))
        print(f"  {body_site}: {n_units:,} genomes")
        for rep in range(n_reps):
            rng = np.random.default_rng(rep)
            perm = labels[rng.permutation(n_units)]
            seen: set = set()
            prev = 0
            for n in subset_sizes:
                seen.update(perm[prev:n].tolist())
                prev = n
                results.append(
                    {
                        "body_site": body_site,
                        "subset": int(n),
                        "replicate": rep,
                        "n_clusters": len(seen),
                    }
                )
    return (
        pl.from_dicts(results)
        .group_by(["body_site", "subset"])
        .agg(
            [
                pl.col("n_clusters").mean().alias("mean_species"),
                pl.col("n_clusters").min().alias("min_species"),
                pl.col("n_clusters").max().alias("max_species"),
            ]
        )
        .sort(["body_site", "subset"])
    )


def _fmt_thousands(x: float, _pos: int) -> str:
    if x == 0:
        return "0"
    if abs(x) >= 1000:
        return f"{x / 1000:.0f}k"
    return f"{x:.0f}"


def plot_by_body_site(
    curve: pl.DataFrame,
    *,
    ylabel: str,
    outfile: str,
    inset_xlim: tuple[int, int] | None = None,
    inset_ylim: tuple[int, int] | None = None,
    legend_lower: bool = False,
) -> None:
    from matplotlib.ticker import FuncFormatter, MaxNLocator

    plt.rcParams.update({"font.size": 18})
    site_order = [s for s in SITE_ORDER if s in curve["body_site"].unique().to_list()]
    small = curve.filter(pl.col("body_site").is_in(["Skin", "Urogenital", "Blood"]))
    if inset_xlim is None:
        inset_xlim = (0, int(small["subset"].max() * 1.08))
    if inset_ylim is None:
        inset_ylim = (0, int(small["max_species"].max() * 1.12))

    fig, ax = plt.subplots(figsize=(6.5, 6))
    for body_site in site_order:
        d = curve.filter(pl.col("body_site") == body_site).sort("subset")
        x = d["subset"].to_numpy()
        y_mean = d["mean_species"].to_numpy()
        y_min = d["min_species"].to_numpy()
        y_max = d["max_species"].to_numpy()
        c = SITE_COLORS[body_site]
        ax.fill_between(x, y_min, y_max, color=c, alpha=0.15, linewidth=0)
        ax.plot(x, y_mean, color=c, lw=2.2, label=body_site)

    ax.set_xlabel("Number of unique genomes sampled", fontdict={"fontweight": "bold"})
    ax.set_ylabel(ylabel, fontdict={"fontweight": "bold"})
    ax.xaxis.set_label_coords(0.42, -0.12)
    ax.grid(True, alpha=0.25)
    ax.set_xlim(0, max(curve["subset"].to_list()) * 1.05)
    ax.set_ylim(0, max(curve["max_species"].to_list()) * 1.05)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.xaxis.set_major_formatter(FuncFormatter(_fmt_thousands))
    ax.yaxis.set_major_formatter(FuncFormatter(_fmt_thousands))
    legend_kwargs = {"title": "Body site", "fontsize": 11, "title_fontsize": 12}
    if legend_lower:
        ax.legend(loc="lower left", **legend_kwargs)
    else:
        ax.legend(loc="upper left", **legend_kwargs)

    axins = ax.inset_axes([0.52, 0.06, 0.42, 0.42])
    for body_site in site_order:
        d = curve.filter(pl.col("body_site") == body_site).sort("subset")
        axins.plot(
            d["subset"].to_numpy(),
            d["mean_species"].to_numpy(),
            color=SITE_COLORS[body_site],
            lw=1.8,
        )
    axins.set_xlim(*inset_xlim)
    axins.set_ylim(*inset_ylim)
    axins.grid(True, alpha=0.2)
    axins.tick_params(labelsize=10)
    axins.xaxis.set_major_formatter(FuncFormatter(_fmt_thousands))
    axins.yaxis.set_major_formatter(FuncFormatter(_fmt_thousands))
    try:
        indicator = ax.indicate_inset_zoom(axins, edgecolor="black", linewidth=2.2)
        if hasattr(indicator, "rectangle"):
            indicator.rectangle.set_linewidth(2.2)
            indicator.rectangle.set_edgecolor("black")
    except Exception:
        pass

    plt.tight_layout()
    out = OUT / outfile
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"wrote {out} (inset x≤{inset_xlim[1]}, y≤{inset_ylim[1]})")


def main() -> None:
    print("Loading unique genomes...")
    genomes = load_unique_genomes()
    species = load_species_table(genomes)
    print(
        f"Species: {species.height:,}; genera: {species['genus_cluster_id'].n_unique():,}; "
        f"families: {species['family_cluster_id'].n_unique():,}; "
        f"TR species: {species.filter(pl.col('is_tr')).height:,}"
    )

    print("Genus accumulation (sampling species)...")
    genus_curve = calculate_accumulation(species, "genus_cluster_id", n_permutations=50)
    genus_curve.write_csv(OUT / "figure_2b_genus_accumulation_r6.tsv", separator="\t")
    plot_overall(
        genus_curve,
        ylabel="Cumulative genera",
        outfile="figure_2b_genus_accumulation_r6.png",
        legend_loc="upper left",
    )

    print("Family accumulation (sampling species)...")
    family_curve = calculate_accumulation(species, "family_cluster_id", n_permutations=50)
    family_curve.write_csv(OUT / "figure_2c_family_accumulation_r6.tsv", separator="\t")
    plot_overall(
        family_curve,
        ylabel="Cumulative families",
        outfile="figure_2c_family_accumulation_r6.png",
        legend_loc="lower right",
    )

    print("Genus by body site...")
    genus_site = body_site_cluster_accumulation(genomes, "genus_cluster_id")
    genus_site.write_csv(OUT / "figure_2b_genus_by_body_site_r6.tsv", separator="\t")
    plot_by_body_site(
        genus_site,
        ylabel="Cumulative genera",
        outfile="figure_2b_genus_by_body_site_r6.png",
    )

    print("Family by body site...")
    family_site = body_site_cluster_accumulation(genomes, "family_cluster_id")
    family_site.write_csv(OUT / "figure_2c_family_by_body_site_r6.tsv", separator="\t")
    plot_by_body_site(
        family_site,
        ylabel="Cumulative families",
        outfile="figure_2c_family_by_body_site_r6.png",
        legend_lower=True,
    )


if __name__ == "__main__":
    main()
