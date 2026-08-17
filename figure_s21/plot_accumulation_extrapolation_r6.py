#!/usr/bin/env python3
"""Separate body-site accumulation plots with iNEXT (Chao) rarefaction/extrapolation."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from matplotlib.ticker import FuncFormatter, MaxNLocator

BASE = Path(__file__).resolve().parent
META = (
    Path(__file__).resolve().parents[1]
    / "figure_1"
    / "uhvdb_v6_metadata_sra.tsv"
)
OUT = BASE / "plots"
OUT.mkdir(exist_ok=True)
R_SCRIPT = BASE / "run_inext_extrapolation.R"

SITE_ORDER = ["Skin", "Gut", "Urogenital", "Airways", "Blood"]
SITE_COLORS = {
    "Skin": "#1f77b4",
    "Gut": "#ff7f0e",
    "Urogenital": "#2ca02c",
    "Airways": "#d62728",
    "Blood": "#9467bd",
}

# iNEXT guidance: species-richness extrapolation is reliable ~up to 2× reference size.
ENDPOINT_MULT = 2.0
# Keep boots modest: Gut/Airways abundance vectors are large (~minutes each).
NBOOT = 5
KNOTS = 40


def _fmt_thousands(x: float, _pos: int) -> str:
    if x == 0:
        return "0"
    if abs(x) >= 1_000_000:
        return f"{x / 1_000_000:.1f}M".replace(".0M", "M")
    if abs(x) >= 1000:
        return f"{x / 1000:.1f}k".replace(".0k", "k")
    return f"{x:.0f}"


def _rscript() -> str:
    env = os.environ.get("RSCRIPT")
    if env:
        return env
    # Prefer micromamba inext env if present
    candidate = Path(
        "/mmfs1/gscratch/pedslabs_hoffman/carsonjm/micromamba_envs/envs/inext/bin/Rscript"
    )
    if candidate.exists():
        return str(candidate)
    return "Rscript"


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
                "species_rep",
                "genus_cluster_id",
                "family_cluster_id",
            ],
            infer_schema_length=5000,
        )
        .filter(pl.col("seq_name") == pl.col("seqhash_rep"))
        .filter(pl.col("species_cluster_id").is_not_null())
    )


def abundance_by_site(genomes: pl.DataFrame, taxon_col: str) -> pl.DataFrame:
    """Individual-based abundances: genomes per taxon within each body site."""
    return (
        genomes.select(["seqhash_rep", "body_site", taxon_col])
        .filter(
            pl.col("body_site").is_not_null()
            & (pl.col("body_site") != "Other")
            & pl.col("body_site").is_in(SITE_ORDER)
            & pl.col(taxon_col).is_not_null()
        )
        .unique()
        .group_by(["body_site", taxon_col])
        .len()
        .rename({taxon_col: "taxon", "len": "abundance"})
        .sort(["body_site", "taxon"])
    )


def run_inext(abundance: pl.DataFrame, out_prefix: Path) -> tuple[pl.DataFrame, pl.DataFrame]:
    ab_path = out_prefix.with_name(out_prefix.name + "_abundance.tsv")
    abundance.write_csv(ab_path, separator="\t")
    cmd = [
        _rscript(),
        str(R_SCRIPT),
        str(ab_path),
        str(out_prefix),
        str(ENDPOINT_MULT),
        str(NBOOT),
        str(KNOTS),
    ]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    curve = pl.read_csv(str(out_prefix) + "_curve.tsv", separator="\t")
    asy = pl.read_csv(str(out_prefix) + "_asy.tsv", separator="\t")
    return curve, asy


def plot_extrapolation(
    curve: pl.DataFrame,
    asy: pl.DataFrame,
    *,
    ylabel: str,
    outfile: str,
    title: str,
) -> pl.DataFrame:
    plt.rcParams.update({"font.size": 12})
    sites = [s for s in SITE_ORDER if s in curve["body_site"].unique().to_list()]
    asy_map = {r["body_site"]: r for r in asy.to_dicts()}

    fig, axes = plt.subplots(1, len(sites), figsize=(3.1 * len(sites), 4.2), sharey=False)
    if len(sites) == 1:
        axes = [axes]

    summary_rows: list[dict] = []

    for ax, site in zip(axes, sites):
        d = curve.filter(pl.col("body_site") == site).sort("m")
        a = asy_map[site]
        c = SITE_COLORS[site]
        x = np.asarray(d["m"].to_numpy(), dtype=float)
        y = np.asarray(d["qD"].to_numpy(), dtype=float)
        methods = d["Method"].to_list()

        # Confidence ribbon when available
        if "qD.LCL" in d.columns and "qD.UCL" in d.columns:
            y_lo = np.asarray(d["qD.LCL"].to_numpy(), dtype=float)
            y_hi = np.asarray(d["qD.UCL"].to_numpy(), dtype=float)
            if np.isfinite(y_lo).any() and np.isfinite(y_hi).any():
                ax.fill_between(x, y_lo, y_hi, color=c, alpha=0.15, linewidth=0)

        # Split rarefaction vs extrapolation
        is_ext = np.array([str(m).lower().startswith("extrap") for m in methods])
        is_obs = np.array([str(m).lower().startswith("observ") for m in methods])
        is_rare = ~(is_ext | is_obs)

        # Continuous rarefaction + observed (solid), extrapolation (dashed)
        rare_mask = is_rare | is_obs
        if rare_mask.any():
            ax.plot(
                x[rare_mask],
                y[rare_mask],
                color=c,
                lw=2.4,
                label="iNEXT rarefaction",
            )
        if is_obs.any():
            ax.plot(
                x[is_obs],
                y[is_obs],
                color="0.15",
                marker="o",
                ms=5,
                lw=0,
                zorder=5,
                label="Observed",
            )
        if is_ext.any():
            # Connect from observed/reference into extrapolation
            join = np.where(rare_mask)[0]
            ext_idx = np.where(is_ext)[0]
            if join.size and ext_idx.size:
                start = join[-1]
                xs = np.concatenate([[x[start]], x[ext_idx]])
                ys = np.concatenate([[y[start]], y[ext_idx]])
            else:
                xs, ys = x[is_ext], y[is_ext]
            ax.plot(xs, ys, color=c, lw=1.8, ls="--", alpha=0.9, label="iNEXT extrapolation")

        s_asy = float(a["S_asy"]) if a["S_asy"] is not None else float("nan")
        if np.isfinite(s_asy):
            ax.axhline(s_asy, color=c, lw=1.2, ls=":", alpha=0.75, label=r"$S_{\mathrm{Chao}}$")

        n_now = float(a["n_now"])
        s_obs = float(a["S_obs"])
        sc_now = float(a["SC_now"]) if a["SC_now"] is not None else float("nan")
        frac = s_obs / s_asy if s_asy > 0 else float("nan")
        endpoint = float(a["endpoint"])
        s_end = float(a["S_endpoint"]) if a["S_endpoint"] is not None else float("nan")

        summary_rows.append(
            {
                "body_site": site,
                "n_now": n_now,
                "S_obs": s_obs,
                "SC_now": sc_now,
                "f1": a.get("f1"),
                "f2": a.get("f2"),
                "S_asy": s_asy,
                "S_asy_se": a.get("S_asy_se"),
                "S_asy_LCL": a.get("S_asy_LCL"),
                "S_asy_UCL": a.get("S_asy_UCL"),
                "frac_obs_of_asy": frac,
                "endpoint": endpoint,
                "S_endpoint": s_end,
                "SC_endpoint": a.get("SC_endpoint"),
            }
        )

        ax.set_title(site, fontweight="bold", fontsize=13)
        ax.grid(True, alpha=0.25)
        ax.set_xlim(0, endpoint * 1.02)
        ymax = np.nanmax([np.nanmax(y), s_asy if np.isfinite(s_asy) else np.nan])
        ax.set_ylim(0, ymax * 1.12 if np.isfinite(ymax) and ymax > 0 else 1.0)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        ax.xaxis.set_major_formatter(FuncFormatter(_fmt_thousands))
        ax.yaxis.set_major_formatter(FuncFormatter(_fmt_thousands))
        ax.tick_params(labelsize=10)

        sc_txt = f"{100.0 * sc_now:.1f}%" if np.isfinite(sc_now) else "NA"
        pct = f"{100.0 * frac:.0f}%" if np.isfinite(frac) else "NA"
        ax.text(
            0.04,
            0.97,
            (
                f"$S_{{\\mathrm{{Chao}}}}$={_fmt_thousands(s_asy, 0)}\n"
                f"obs={pct} of Chao\n"
                f"SC={sc_txt}"
            ),
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=9,
            linespacing=1.35,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="0.8", alpha=0.9),
        )

    axes[0].set_ylabel(ylabel, fontdict={"fontweight": "bold"})
    fig.supxlabel("Number of unique genomes sampled", fontweight="bold", y=0.02)

    # Shared legend from first panel that has all handle types
    wanted = [
        "Observed",
        "iNEXT rarefaction",
        "iNEXT extrapolation",
        r"$S_{\mathrm{Chao}}$",
    ]
    handles, labels = [], []
    for ax in axes:
        h, l = ax.get_legend_handles_labels()
        for hh, ll in zip(h, l):
            if ll in wanted and ll not in labels:
                handles.append(hh)
                labels.append(ll)
    order = [labels.index(l) for l in wanted if l in labels]
    fig.legend(
        [handles[i] for i in order],
        [labels[i] for i in order],
        loc="upper center",
        ncol=4,
        frameon=False,
        fontsize=10,
        bbox_to_anchor=(0.5, 1.02),
    )
    fig.suptitle(title, fontweight="bold", y=1.08, fontsize=14)
    fig.tight_layout()
    out = OUT / outfile
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")

    fit_df = pl.DataFrame(summary_rows)
    tsv = out.with_suffix(".tsv")
    fit_df.write_csv(tsv, separator="\t")
    print(f"wrote {tsv}")
    return fit_df


def main() -> None:
    print("Loading unique genomes...")
    genomes = load_unique_genomes()

    jobs = [
        (
            "species_rep",
            "Cumulative species",
            "figure_2a_species_by_body_site_extrapolation_r6",
            "Species accumulation — iNEXT rarefaction / extrapolation",
        ),
        (
            "genus_cluster_id",
            "Cumulative genera",
            "figure_2b_genus_by_body_site_extrapolation_r6",
            "Genus accumulation — iNEXT rarefaction / extrapolation",
        ),
        (
            "family_cluster_id",
            "Cumulative families",
            "figure_2c_family_by_body_site_extrapolation_r6",
            "Family accumulation — iNEXT rarefaction / extrapolation",
        ),
    ]

    for taxon_col, ylabel, stem, title in jobs:
        print(f"\n=== {stem} ({taxon_col}) ===")
        ab = abundance_by_site(genomes, taxon_col)
        print(
            ab.group_by("body_site")
            .agg(
                [
                    pl.col("abundance").sum().alias("n"),
                    pl.len().alias("S"),
                ]
            )
            .sort("body_site")
        )
        out_prefix = OUT / stem
        curve, asy = run_inext(ab, out_prefix)
        plot_extrapolation(
            curve,
            asy,
            ylabel=ylabel,
            outfile=f"{stem}.png",
            title=title,
        )


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as e:
        print(f"iNEXT Rscript failed with exit {e.returncode}", file=sys.stderr)
        raise
