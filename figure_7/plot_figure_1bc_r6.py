#!/usr/bin/env python3
"""Regenerate figure 1b (release sizes) and 1c (pies) for UHVDB r6."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import seaborn as sns

BASE = Path(__file__).resolve().parent
MANUSCRIPT_FIG1 = Path(
    "/mmfs1/gscratch/pedslabs_hoffman/carsonjm/CFPhageome/repos/UHVDB/"
    "uhvdb-manuscript/figure_1"
)
OUT = BASE / "plots"
OUT.mkdir(exist_ok=True)

META = BASE / "uhvdb_v6_metadata_sra.tsv"

RELEASE_SEQHASHERS = {
    "r1": MANUSCRIPT_FIG1 / "uhgv_hq_hc_results/uhvdb_2026-03-23/uhvdb_seqhasher.tsv.gz",
    "r2": MANUSCRIPT_FIG1 / "uhgv_hq_lc_results/uhvdb_2026-03-25/uhvdb_seqhasher.tsv.gz",
    "r3": MANUSCRIPT_FIG1 / "uhvdb_virus_db_results/uhvdb_2026-03-26/uhvdb_seqhasher.tsv.gz",
    "r4": MANUSCRIPT_FIG1 / "uhvdb_human_metag_results/uhvdb_2026-03-26-2/uhvdb_seqhasher.tsv.gz",
}


def plot_1b_release_sizes() -> None:
    """Unique sequence hashes after each release (r1–r4 from seqhashers; r5/r6 from metadata)."""
    counts: dict[str, int] = {}
    for release, path in RELEASE_SEQHASHERS.items():
        n = pl.read_csv(path, separator="\t", columns=["hash"]).n_unique("hash")
        counts[release] = n
        print(f"{release}: {n:,} unique sequences")

    meta = pl.scan_csv(META, separator="\t").select(["hash", "added_in_release"]).collect()
    # Cumulative unique hashes through r5 and r6 from metadata release labels
    seen: set[str] = set()
    for release in ["r1", "r2", "r3", "r4", "r5", "r6"]:
        seen |= set(
            meta.filter(pl.col("added_in_release") == release)["hash"]
            .drop_nulls()
            .to_list()
        )
        if release in ("r5", "r6"):
            counts[release] = len(seen)
            print(f"{release}: {len(seen):,} unique sequences (cumulative from metadata)")

    plot_releases = ["r4", "r5", "r6"]
    plot_df = pl.DataFrame(
        {
            "Release": plot_releases,
            "Sequence count": [counts[r] for r in plot_releases],
        }
    )

    plt.rcParams.update({"font.size": 14})
    plt.figure(figsize=(5, 5))
    sns.barplot(
        data=plot_df.to_pandas(),
        x="Release",
        y="Sequence count",
        color="#2E86AB",
    )
    plt.ylabel("Unique sequences", fontdict={"fontweight": "bold"})
    plt.xlabel("UHVDB release", fontdict={"fontweight": "bold"})
    for i, row in enumerate(plot_df.iter_rows(named=True)):
        plt.text(
            i,
            row["Sequence count"],
            f"{row['Sequence count']:,}",
            ha="center",
            va="bottom",
            fontsize=11,
        )
    plt.tight_layout()
    out = OUT / "figure_1b_release_sizes_r6.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"wrote {out}")


def load_unique_viruses() -> pl.DataFrame:
    return (
        pl.read_csv(
            META,
            separator="\t",
            columns=[
                "seq_name",
                "source_db",
                "body_site",
                "db_type",
                "checkv_quality",
                "seqhash_rep",
            ],
            infer_schema_length=5000,
        )
        .filter(~pl.col("source_db").is_in(["MOTUS_DB", "CF_METAG"]))
        .filter(pl.col("seq_name") == pl.col("seqhash_rep"))
    )


def _save_pie(fig: plt.Figure, name: str) -> None:
    out = OUT / name
    fig.savefig(out, dpi=300, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"wrote {out}")


def plot_1c_source_db(df: pl.DataFrame) -> None:
    source_db_plot = (
        df.with_columns(
            pl.when(
                pl.col("source_db").is_in(
                    ["OVD", "SMGC", "VMGC", "LOGAN", "CHVD", "CNGVR", "METAVR"]
                )
            )
            .then(pl.lit("Other"))
            .otherwise(pl.col("source_db"))
            .alias("source_db")
        )
        .group_by("source_db")
        .len()
        .sort("len", descending=True)
    )
    print("source_db:", source_db_plot)

    plt.rcParams.update({"font.size": 20})
    fig, ax = plt.subplots(figsize=(7, 7))
    fig.patch.set_alpha(0)
    ax.set_facecolor("none")
    wedges, _ = ax.pie(source_db_plot["len"], startangle=90)

    r0, r_step, min_dy = 1.15, 0.08, 0.22
    pts = []
    for w, label in zip(wedges, source_db_plot["source_db"].to_list()):
        th = np.deg2rad((w.theta1 + w.theta2) / 2)
        pts.append(
            {
                "label": label,
                "th": th,
                "x_rim": float(np.cos(th)),
                "y_rim": float(np.sin(th)),
                "side": "right" if np.cos(th) >= 0 else "left",
                "r": r0,
            }
        )
    for side in ("right", "left"):
        side_pts = sorted([p for p in pts if p["side"] == side], key=lambda p: p["th"])
        for _ in range(100):
            for p in side_pts:
                p["x"] = p["r"] * p["x_rim"]
                p["y"] = p["r"] * p["y_rim"]
            moved = False
            for i in range(1, len(side_pts)):
                a, b = side_pts[i - 1], side_pts[i]
                if abs(a["y"] - b["y"]) < min_dy and abs(a["x"] - b["x"]) < 0.9:
                    b["r"] += r_step
                    moved = True
            if not moved:
                break
    for p in pts:
        ax.plot([p["x_rim"], p["x"]], [p["y_rim"], p["y"]], color="black", lw=0.8, clip_on=False)
        ax.text(
            p["x"] + (0.04 if p["side"] == "right" else -0.04),
            p["y"],
            p["label"],
            ha="left" if p["side"] == "right" else "right",
            va="center",
            fontsize=20,
            clip_on=False,
        )
    ax.set_xlim(-1.9, 1.9)
    ax.set_ylim(-1.9, 1.9)
    ax.set_aspect("equal")
    ax.axis("off")
    plt.tight_layout()
    _save_pie(fig, "figure_1c_source_db_r6.png")


def plot_1c_db_type(df: pl.DataFrame) -> None:
    source_type_plot = df.group_by("db_type").len().sort("len", descending=True)
    print("db_type:", source_type_plot)
    plt.rcParams.update({"font.size": 24})
    fig, ax = plt.subplots(figsize=(7, 7))
    fig.patch.set_alpha(0)
    ax.set_facecolor("none")
    wedges, _ = ax.pie(source_type_plot["len"], startangle=90)
    r_stub = 1.12
    for w, label in zip(wedges, source_type_plot["db_type"].to_list()):
        th = np.deg2rad((w.theta1 + w.theta2) / 2)
        x_rim, y_rim = np.cos(th), np.sin(th)
        x, y = r_stub * x_rim, r_stub * y_rim
        side_right = x_rim >= 0
        ax.plot([x_rim, x], [y_rim, y], color="black", lw=0.8, clip_on=False)
        ax.text(
            x + (0.04 if side_right else -0.04),
            y,
            label,
            ha="left" if side_right else "right",
            va="center",
            fontsize=24,
            clip_on=False,
        )
    ax.set_xlim(-1.9, 1.9)
    ax.set_ylim(-1.9, 1.9)
    ax.set_aspect("equal")
    ax.axis("off")
    plt.tight_layout()
    _save_pie(fig, "figure_1c_db_type_r6.png")


def plot_1c_body_site(df: pl.DataFrame) -> None:
    body_site_plot = df.group_by("body_site").len().sort("len", descending=True)
    print("body_site:", body_site_plot)
    plt.rcParams.update({"font.size": 24})
    fig, ax = plt.subplots(figsize=(7, 7))
    fig.patch.set_alpha(0)
    ax.set_facecolor("none")
    wedges, _ = ax.pie(body_site_plot["len"], startangle=90)
    r_line, x_pad = 1.06, 0.003
    label_radius_add = {"Urogenital": 0.25, "Skin": 0.1, "Other": 0.1, "Blood": 0.15}
    pts = []
    for w, label in zip(wedges, body_site_plot["body_site"].to_list()):
        th = np.deg2rad((w.theta1 + w.theta2) / 2)
        x, y = np.cos(th), np.sin(th)
        pts.append(
            {
                "label": label,
                "x": x,
                "y": y,
                "r_end": r_line + label_radius_add.get(label, 0.0),
                "side": "right" if x >= 0 else "left",
            }
        )
    for p in pts:
        x_end = p["r_end"] * p["x"]
        y_end = p["r_end"] * p["y"]
        ax.plot([p["x"], x_end], [p["y"], y_end], color="black", lw=0.8)
        x_txt = x_end + x_pad if p["side"] == "right" else x_end - x_pad
        ax.text(x_txt, y_end, p["label"], ha="left" if p["side"] == "right" else "right", va="center", fontsize=24)
    ax.set_xlim(-1.9, 1.9)
    ax.set_ylim(-1.9, 1.9)
    ax.set_aspect("equal")
    plt.tight_layout()
    _save_pie(fig, "figure_1c_body_site_r6.png")


def plot_1c_checkv(df: pl.DataFrame) -> None:
    checkv_quality_plot = df.group_by("checkv_quality").len().sort("len", descending=True)
    print("checkv_quality:", checkv_quality_plot)
    plt.rcParams.update({"font.size": 24})
    fig, ax = plt.subplots(figsize=(7, 7))
    fig.patch.set_alpha(0)
    ax.set_facecolor("none")
    wedges, _ = ax.pie(checkv_quality_plot["len"], startangle=90)
    v_len, x_pad = 0.14, 0.01
    pts = []
    for w, label in zip(wedges, checkv_quality_plot["checkv_quality"].to_list()):
        th = np.deg2rad((w.theta1 + w.theta2) / 2)
        x, y = np.cos(th), np.sin(th)
        pts.append(
            {
                "label": label,
                "x": x,
                "y": y,
                "y_anchor": y + (v_len if y >= 0 else -v_len),
                "side": "right" if x >= 0 else "left",
            }
        )
    for p in pts:
        ax.plot([p["x"], p["x"]], [p["y"], p["y_anchor"]], color="black", lw=0.8)
        x_txt = p["x"] + x_pad if p["side"] == "right" else p["x"] - x_pad
        ax.text(x_txt, p["y_anchor"], p["label"], ha="left" if p["side"] == "right" else "right", va="center", fontsize=24)
    ax.set_xlim(-1.9, 1.9)
    ax.set_ylim(-1.9, 1.9)
    ax.set_aspect("equal")
    plt.tight_layout()
    _save_pie(fig, "figure_1c_checkv_quality_r6.png")


def main() -> None:
    print("=== Figure 1b ===")
    plot_1b_release_sizes()
    print("=== Figure 1c ===")
    df = load_unique_viruses()
    print(f"Unique viruses (seq_name == seqhash_rep): {df.height:,}")
    plot_1c_source_db(df)
    plot_1c_db_type(df)
    plot_1c_body_site(df)
    plot_1c_checkv(df)


if __name__ == "__main__":
    main()
