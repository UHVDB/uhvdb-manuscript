#!/usr/bin/env python3
"""Regenerate figure 2a species accumulation curves for UHVDB r6."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import seaborn as sns

BASE = Path(__file__).resolve().parent
META = (
    Path(__file__).resolve().parents[1]
    / "figure_1"
    / "uhvdb_v6_metadata_sra.tsv"
)
OUT = BASE / "plots"
OUT.mkdir(exist_ok=True)


def load_unique_genomes() -> pl.DataFrame:
    """One row per unique virus (seqhash representative)."""
    return (
        pl.read_csv(
            META,
            separator="\t",
            columns=[
                "seq_name",
                "seqhash_rep",
                "body_site",
                "species_cluster_id",
                "species_rep",
                "completeness_method",
            ],
            infer_schema_length=5000,
        )
        .filter(pl.col("seq_name") == pl.col("seqhash_rep"))
        .filter(pl.col("species_cluster_id").is_not_null())
        .with_columns(pl.col("completeness_method").fill_null(""))
    )


def accumulate_once(
    species: np.ndarray, is_tr: np.ndarray, rng: np.random.Generator
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (total_species, non_singletons, tr_species) cumulative curves."""
    import pandas as pd

    order = rng.permutation(species.shape[0])
    s = pd.Series(species[order])
    t = is_tr[order]

    total = (~s.duplicated()).cumsum().to_numpy(dtype=np.int32)
    non_sing = (s.groupby(s, sort=False).cumcount() == 1).cumsum().to_numpy(dtype=np.int32)
    s_tr = s.where(t)
    tr_sp = (s_tr.notna() & ~s_tr.duplicated()).cumsum().to_numpy(dtype=np.int32)
    return total, non_sing, tr_sp
def calculate_complex_accumulation(
    df: pl.DataFrame, n_permutations: int = 50
) -> pl.DataFrame:
    species = df["species_cluster_id"].to_numpy()
    is_tr = df["completeness_method"].str.contains("TR").fill_null(False).to_numpy()
    n = species.shape[0]
    sum_total = np.zeros(n, dtype=np.float64)
    sum_non = np.zeros(n, dtype=np.float64)
    sum_tr = np.zeros(n, dtype=np.float64)

    for i in range(n_permutations):
        rng = np.random.default_rng(i)
        total, non_sing, tr_sp = accumulate_once(species, is_tr, rng)
        sum_total += total
        sum_non += non_sing
        sum_tr += tr_sp
        if (i + 1) % 10 == 0:
            print(f"  permutation {i + 1}/{n_permutations}")

    return pl.DataFrame(
        {
            "genomes_sampled": np.arange(1, n + 1),
            "mean_total": sum_total / n_permutations,
            "mean_non_singletons": sum_non / n_permutations,
            "mean_tr": sum_tr / n_permutations,
        }
    )


def plot_overall(accumulation_curve: pl.DataFrame) -> None:
    sns.reset_orig()
    plt.rcParams.update({"font.size": 14})
    # Thin for plotting speed (curve is smooth)
    step = max(1, accumulation_curve.height // 5000)
    plot_df = accumulation_curve.gather_every(step)
    x = plot_df["genomes_sampled"].to_numpy()
    y_total = plot_df["mean_total"].to_numpy()
    y_non_sing = plot_df["mean_non_singletons"].to_numpy()

    plt.figure(figsize=(7, 6))
    plt.grid(True, alpha=0.3)
    sns.lineplot(x=x, y=y_total, color="#2E86AB", linewidth=2.5, label="Total")
    sns.lineplot(x=x, y=y_non_sing, color="#D62828", linewidth=2.5, label="Non-singleton")
    plt.xlabel("Number of genomes sampled", fontdict={"fontweight": "bold"})
    plt.ylabel("Cumulative species", fontdict={"fontweight": "bold"})
    plt.xlim(0, max(x))
    plt.ylim(0, max(y_total) * 1.05)
    plt.legend(title="Cluster type", loc="upper left")
    plt.tight_layout()
    out = OUT / "figure_2a_species_accumulation_r6.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"wrote {out}")


def body_site_accumulation(
    df: pl.DataFrame, n_reps: int = 50, n_points: int = 100
) -> pl.DataFrame:
    body_site_species_df = (
        df.select(["seqhash_rep", "body_site", "species_rep"])
        .filter(
            pl.col("body_site").is_not_null()
            & (pl.col("body_site") != "Other")
            & pl.col("species_rep").is_not_null()
        )
        .unique()
    )
    results = []
    for body_site in body_site_species_df["body_site"].unique().to_list():
        site_df = body_site_species_df.filter(pl.col("body_site") == body_site)
        n_units = site_df.height
        if n_units == 0:
            continue
        species = site_df["species_rep"].to_numpy()
        subset_sizes = np.unique(np.linspace(1, n_units, min(n_points, n_units), dtype=int))
        print(f"  {body_site}: {n_units:,} genomes")
        for rep in range(n_reps):
            rng = np.random.default_rng(rep)
            perm = species[rng.permutation(n_units)]
            seen = set()
            prev = 0
            for n in subset_sizes:
                seen.update(perm[prev:n].tolist())
                prev = n
                results.append(
                    {
                        "body_site": body_site,
                        "subset": int(n),
                        "replicate": rep,
                        "species": len(seen),
                    }
                )
    return (
        pl.from_dicts(results)
        .group_by(["body_site", "subset"])
        .agg(
            [
                pl.col("species").mean().alias("mean_species"),
                pl.col("species").min().alias("min_species"),
                pl.col("species").max().alias("max_species"),
            ]
        )
        .sort(["body_site", "subset"])
    )


def _fmt_thousands(x: float, _pos: int) -> str:
    """Compact tick labels: 0, 100k, 200k, …"""
    if x == 0:
        return "0"
    if abs(x) >= 1000:
        return f"{x / 1000:.0f}k"
    return f"{x:.0f}"


def plot_by_body_site(curve: pl.DataFrame) -> None:
    from matplotlib.ticker import FuncFormatter, MaxNLocator

    plt.rcParams.update({"font.size": 18})
    site_order = [
        s
        for s in ["Skin", "Gut", "Urogenital", "Airways", "Blood"]
        if s in curve["body_site"].unique().to_list()
    ]
    colors = {
        "Skin": "#1f77b4",
        "Gut": "#ff7f0e",
        "Urogenital": "#2ca02c",
        "Airways": "#d62728",
        "Blood": "#9467bd",
    }
    # Inset must cover full Skin / Urogenital / Blood curves
    small_sites = ["Skin", "Urogenital", "Blood"]
    small = curve.filter(pl.col("body_site").is_in(small_sites))
    inset_xmax = int(small["subset"].max() * 1.08)
    inset_ymax = int(small["max_species"].max() * 1.12)

    fig, ax = plt.subplots(figsize=(6.5, 6))
    for body_site in site_order:
        d = curve.filter(pl.col("body_site") == body_site).sort("subset")
        x = d["subset"].to_numpy()
        y_mean = d["mean_species"].to_numpy()
        y_min = d["min_species"].to_numpy()
        y_max = d["max_species"].to_numpy()
        c = colors[body_site]
        ax.fill_between(x, y_min, y_max, color=c, alpha=0.15, linewidth=0)
        ax.plot(x, y_mean, color=c, lw=2.2, label=body_site)

    ax.set_xlabel("Number of unique genomes sampled", fontdict={"fontweight": "bold"})
    ax.set_ylabel("Cumulative species", fontdict={"fontweight": "bold"})
    ax.xaxis.set_label_coords(0.42, -0.12)
    ax.legend(title="Body site", loc="upper left", fontsize=11, title_fontsize=12)
    ax.grid(True, alpha=0.25)
    ax.set_xlim(0, max(curve["subset"].to_list()) * 1.05)
    ax.set_ylim(0, max(curve["max_species"].to_list()) * 1.05)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.xaxis.set_major_formatter(FuncFormatter(_fmt_thousands))
    ax.yaxis.set_major_formatter(FuncFormatter(_fmt_thousands))

    axins = ax.inset_axes([0.52, 0.06, 0.42, 0.42])
    for body_site in site_order:
        d = curve.filter(pl.col("body_site") == body_site).sort("subset")
        axins.plot(d["subset"].to_numpy(), d["mean_species"].to_numpy(), color=colors[body_site], lw=1.8)
    axins.set_xlim(0, inset_xmax)
    axins.set_ylim(0, inset_ymax)
    axins.grid(True, alpha=0.2)
    axins.tick_params(labelsize=10)
    axins.xaxis.set_major_formatter(FuncFormatter(_fmt_thousands))
    axins.yaxis.set_major_formatter(FuncFormatter(_fmt_thousands))
    rect, connectors = ax.indicate_inset_zoom(axins, edgecolor="black", linewidth=2.5)
    rect.set_linewidth(2.5)
    for conn in connectors:
        conn.set_linewidth(2.0)
        conn.set_color("black")

    plt.tight_layout()
    out = OUT / "figure_2a_species_by_body_site_r6.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"wrote {out} (inset x≤{inset_xmax}, y≤{inset_ymax})")


def main() -> None:
    print("Loading unique genomes...")
    df = load_unique_genomes()
    print(
        f"Unique genomes: {df.height:,}; species: {df['species_cluster_id'].n_unique():,}; "
        f"TR genomes: {df.filter(pl.col('completeness_method').str.contains('TR')).height:,}"
    )

    print("Overall species accumulation (50 permutations)...")
    curve = calculate_complex_accumulation(df, n_permutations=50)
    curve.write_csv(OUT / "figure_2a_species_accumulation_r6.tsv", separator="\t")
    plot_overall(curve)

    print("Body-site species accumulation...")
    site_curve = body_site_accumulation(df, n_reps=50, n_points=100)
    site_curve.write_csv(OUT / "figure_2a_species_by_body_site_r6.tsv", separator="\t")
    plot_by_body_site(site_curve)


if __name__ == "__main__":
    main()
