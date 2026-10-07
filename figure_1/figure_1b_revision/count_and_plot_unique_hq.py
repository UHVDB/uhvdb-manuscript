#!/usr/bin/env python3
"""Count unique HQ sequence hashes and plot stacked human vs non-human bars."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import polars as pl
from matplotlib.patches import Patch

BASE = Path("/mmfs1/gscratch/pedslabs_hoffman/carsonjm/CFPhageome/repos/UHVDB")
WORKDIR = BASE / "uhvdb-manuscript/figure_1/figure_1b_revision"
DATA = WORKDIR / "data"
PLOTS = WORKDIR / "plots"
PLOTS.mkdir(parents=True, exist_ok=True)

UHVDB_META = BASE / "toolkit2/databases/uhvdb/v6/analyze/uhvdb_metadata.tsv.gz"


def load_seqhasher_dir(chunk_dir: Path) -> pl.DataFrame:
    """Load id + hash from all chunk TSVs (no header: id, hash, sequence)."""
    files = sorted(chunk_dir.glob("*.seqhasher.tsv.gz"))
    if not files:
        raise FileNotFoundError(f"No seqhasher chunks in {chunk_dir}")
    frames = []
    for f in files:
        df = pl.read_csv(
            f,
            separator="\t",
            has_header=False,
            columns=[0, 1],
            new_columns=["seq_id", "hash"],
            infer_schema_length=0,
        )
        frames.append(df)
        print(f"  {f.name}: {df.height:,} rows")
    out = pl.concat(frames)
    print(f"  total rows: {out.height:,}; unique hashes: {out.n_unique('hash'):,}")
    return out


def main() -> None:
    counts: dict[str, dict[str, int]] = {}

    # --- UHGV: all human-associated ---
    print("=== UHGV ===")
    uhgv = load_seqhasher_dir(DATA / "seqhasher_chunks/uhgv")
    n = uhgv.n_unique("hash")
    counts["UHGV"] = {"total": n, "human": n}

    # --- metaVR ---
    print("=== metaVR ===")
    metavr = load_seqhasher_dir(DATA / "seqhasher_chunks/metavr").with_columns(
        pl.col("seq_id").str.split("|").list.first().alias("seq_id")
    )
    metavr_human = set(
        pl.read_csv(
            DATA / "metavr_hq_human_ids.txt",
            has_header=False,
            new_columns=["id"],
        )["id"].to_list()
    )
    total = metavr.n_unique("hash")
    human = metavr.filter(pl.col("seq_id").is_in(metavr_human)).n_unique("hash")
    counts["metaVR"] = {"total": total, "human": human}
    print(f"  human-associated unique hashes: {human:,}")

    # --- VIRE ---
    print("=== VIRE ===")
    vire = load_seqhasher_dir(DATA / "seqhasher_chunks/vire")
    vire_human = set(
        pl.read_csv(
            DATA / "vire_hq_human_ids.txt",
            has_header=False,
            new_columns=["id"],
        )["id"].to_list()
    )
    total = vire.n_unique("hash")
    human = vire.filter(pl.col("seq_id").is_in(vire_human)).n_unique("hash")
    counts["VIRE"] = {"total": total, "human": human}
    print(f"  human-associated unique hashes: {human:,}")

    # --- UHVDB r6: already HQ; all human-associated ---
    print("=== UHVDB r6 ===")
    uhvdb = (
        pl.scan_csv(UHVDB_META, separator="\t")
        .filter(pl.col("checkv_quality").is_in(["Complete", "High-quality"]))
        .select("hash")
        .unique()
        .collect()
    )
    n = uhvdb.height
    counts["UHVDB"] = {"total": n, "human": n}
    print(f"  unique hashes: {n:,}")

    order = ["UHGV", "metaVR", "VIRE", "UHVDB"]
    plot_df = pl.DataFrame(
        {
            "database": order,
            "total": [counts[k]["total"] for k in order],
            "human": [counts[k]["human"] for k in order],
        }
    ).with_columns((pl.col("total") - pl.col("human")).alias("non_human"))
    print(plot_df)

    out_tsv = DATA / "unique_hq_hash_counts.tsv"
    plot_df.write_csv(out_tsv, separator="\t")
    print(f"wrote {out_tsv}")
    plot_unique_hq_bars(plot_df)


def plot_unique_hq_bars(plot_df: pl.DataFrame | None = None) -> Path:
    """Stacked bars: per-database hue; human = dark + diagonal hatch; other = light."""
    if plot_df is None:
        plot_df = pl.read_csv(DATA / "unique_hq_hash_counts.tsv", separator="\t")

    # Distinct hue per database (dark = human, light = other)
    colors = {
        "UHGV": {"human": "#1B4F72", "other": "#AED6F1"},
        "metaVR": {"human": "#196F3D", "other": "#ABEBC6"},
        "VIRE": {"human": "#922B21", "other": "#F5B7B1"},
        "UHVDB": {"human": "#6C3483", "other": "#D2B4DE"},
    }

    databases = plot_df["database"].to_list()
    human = plot_df["human"].to_list()
    other = plot_df["non_human"].to_list()
    totals = plot_df["total"].to_list()
    x = list(range(len(databases)))

    plt.rcParams.update({"font.size": 14})
    fig, ax = plt.subplots(figsize=(6, 7))

    # Human-associated (bottom): dark + diagonal hatch
    for i, db in enumerate(databases):
        ax.bar(
            x[i],
            human[i],
            width=0.7,
            color=colors[db]["human"],
            hatch="///",
            edgecolor="white",
            linewidth=0.6,
            label="From human metagenome" if i == 0 else None,
        )
    # Other (top): light
    for i, db in enumerate(databases):
        ax.bar(
            x[i],
            other[i],
            width=0.7,
            bottom=human[i],
            color=colors[db]["other"],
            edgecolor="white",
            linewidth=0.6,
            label="Other" if i == 0 else None,
        )

    ax.set_ylabel("Unique, high-quality viruses", fontdict={"fontweight": "bold"})
    ax.set_xlabel("Database", fontdict={"fontweight": "bold"})
    ax.set_xticks(x)
    ax.set_xticklabels(databases)

    ymax = max(totals)
    tick_max = int(((ymax // 200_000) + 1) * 200_000)
    ticks = list(range(0, tick_max + 1, 200_000))
    ax.set_yticks(ticks)
    ax.set_ylim(0, tick_max)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{int(v):,}"))

    # Legend: hatch/style for human vs other (neutral gray so it isn't tied to one DB)
    legend_handles = [
        Patch(facecolor="#555555", hatch="///", edgecolor="white", label="From human metagenome"),
        Patch(facecolor="#CCCCCC", edgecolor="white", label="Other"),
    ]
    ax.legend(handles=legend_handles, frameon=False, loc="upper left")

    for i, total in enumerate(totals):
        ax.text(i, total, f"{total:,}", ha="center", va="bottom", fontsize=11)

    fig.tight_layout()
    out_png = PLOTS / "figure_1b_unique_hq_by_database.png"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_png}")
    return out_png


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == "--plot-only":
        plot_unique_hq_bars()
    else:
        main()
