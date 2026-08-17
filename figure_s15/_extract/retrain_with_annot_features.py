#!/usr/bin/env python3
"""Retrain figure_s15 Liang inactive-virus RF with added PHROG/Bakta annot features."""
from __future__ import annotations

import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline

FIG = Path(__file__).resolve().parents[1]
ROOT = FIG.parents[1]
TRAIN = FIG / "liang_train_enrich_v_unenrich_final.tsv"
META = ROOT / "uhvdb-manuscript/figure_s18/uhvdb_v5_final_metadata_v2.tsv.gz"
ANN = ROOT / "uhvdb-manuscript/figure_s18/uhvdb_v5_protein_annotations.tsv.gz"
PHAROKKA = (
    ROOT
    / "uhvdb-manuscript/figure_s18/uhvdb_cf_results/uhvdb_2026-04-03/uhvdb_pharokka.tsv.gz"
)
OPT = FIG / "optimal_features.tsv"
OUT_DIR = FIG / "annot_feature_ablation"
OUT_DIR.mkdir(parents=True, exist_ok=True)

NEW_FEATURES = [
    "annot_has_methylase",
    "annot_has_PHROG8",
    "annot_has_repressor_excisionase",  # PHROG 6928 or 2790
    "annot_has_IS",
    "annot_has_head_packaging",
    "annot_has_excisionase_regulator",  # pharokka annot string
]


def build_species_annot_features(species_ids: set[int]) -> pd.DataFrame:
    """Binary annotation features on species representatives."""
    meta = pl.read_csv(META, separator="\t", infer_schema_length=10000)
    reps = (
        meta.filter(pl.col("seq_name") == pl.col("seqhash_rep"))
        .with_columns(pl.col("species_cluster_id").cast(pl.Int64))
        .filter(pl.col("species_cluster_id").is_in(list(species_ids)))
        .select(
            [
                "species_cluster_id",
                pl.col("uhvdb_id").alias("genomovar_rep"),
            ]
        )
        .unique("species_cluster_id", keep="first")
    )
    rep_ids = set(reps["genomovar_rep"].to_list())
    sp_of = dict(zip(reps["genomovar_rep"].to_list(), reps["species_cluster_id"].to_list()))
    print(f"species reps to annotate: {len(rep_ids)}")

    hash_to_reps = defaultdict(set)
    feats = defaultdict(set)  # genomovar -> feature names

    print(f"scanning {ANN}")
    with gzip.open(ANN, "rt") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rep = row["genomovar_rep"]
            if rep not in rep_ids:
                continue
            h = (row.get("hash") or "").strip()
            if h:
                hash_to_reps[h].add(rep)
            annot = (row.get("pharokka_annot") or "").strip()
            al = annot.lower()
            if al == "dna methyltransferase" or "methyltransferase" in al or "methylase" in al:
                feats[rep].add("annot_has_methylase")
            if al == "excisionase and transcriptional regulator":
                feats[rep].add("annot_has_excisionase_regulator")
            cat = (row.get("pharokka_category") or "").strip()
            if cat == "head and packaging":
                feats[rep].add("annot_has_head_packaging")
            bakta = (row.get("bakta_acc") or "").strip()
            if bakta.startswith("IS:"):
                feats[rep].add("annot_has_IS")

    needed = set(hash_to_reps)
    print(f"scanning pharokka for PHROG IDs ({len(needed)} hashes)")
    with gzip.open(PHAROKKA, "rt") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["ID"] not in needed:
                continue
            pid = row["phrog"]
            if pid == "8":
                for rep in hash_to_reps[row["ID"]]:
                    feats[rep].add("annot_has_PHROG8")
            elif pid in {"6928", "2790"}:
                for rep in hash_to_reps[row["ID"]]:
                    feats[rep].add("annot_has_repressor_excisionase")
            # also catch methylase via pharokka annot if missed
            annot = (row.get("annot") or "").strip().lower()
            if "methyltransferase" in annot or "methylase" in annot:
                for rep in hash_to_reps[row["ID"]]:
                    feats[rep].add("annot_has_methylase")

    rows = []
    for rep, sp in sp_of.items():
        fset = feats.get(rep, set())
        row = {"species_cluster_id": int(sp)}
        for ft in NEW_FEATURES:
            row[ft] = 1.0 if ft in fset else 0.0
        rows.append(row)
    out = pd.DataFrame(rows)
    print("annot prevalence (species reps):")
    for ft in NEW_FEATURES:
        print(f"  {ft}: {out[ft].sum():.0f}/{len(out)}")
    return out


def run_cv(df: pd.DataFrame, feature_cols: list[str], title: str) -> dict:
    d = df[df["pr_cat"].isin(["TP", "FP"])].copy()
    y = (d["pr_cat"] == "TP").astype(int)
    groups = d["group"]
    X = d.reindex(columns=feature_cols)
    missing = [c for c in feature_cols if c not in d.columns]
    if missing:
        raise KeyError(missing)

    pipe = Pipeline(
        steps=[
            (
                "preprocessor",
                ColumnTransformer(
                    transformers=[("num", SimpleImputer(strategy="median"), feature_cols)]
                ),
            ),
            (
                "classifier",
                RandomForestClassifier(
                    n_estimators=100,
                    class_weight="balanced",
                    random_state=42,
                    n_jobs=-1,
                ),
            ),
        ]
    )
    gkf = GroupKFold(n_splits=min(5, d["group"].nunique()))
    all_y, all_p, all_hat = [], [], []
    fold_rows = []
    importances = np.zeros(len(feature_cols))
    n_imp = 0
    for fold, (tr, te) in enumerate(gkf.split(X, y, groups), start=1):
        pipe.fit(X.iloc[tr], y.iloc[tr])
        proba = pipe.predict_proba(X.iloc[te])[:, 1]
        pred = pipe.predict(X.iloc[te])
        yt = y.iloc[te]
        all_y.extend(yt.tolist())
        all_p.extend(proba.tolist())
        all_hat.extend(pred.tolist())
        auroc = roc_auc_score(yt, proba) if yt.nunique() > 1 else np.nan
        auprc = average_precision_score(yt, proba) if yt.nunique() > 1 else np.nan
        mcc = matthews_corrcoef(yt, pred) if yt.nunique() > 1 else np.nan
        fold_rows.append({"fold": fold, "AUROC": auroc, "AUPRC": auprc, "MCC": mcc})
        clf = pipe.named_steps["classifier"]
        importances += clf.feature_importances_
        n_imp += 1

    pooled = {
        "title": title,
        "n_features": len(feature_cols),
        "n": len(d),
        "n_pos": int(y.sum()),
        "AUROC": roc_auc_score(all_y, all_p),
        "AUPRC": average_precision_score(all_y, all_p),
        "MCC": matthews_corrcoef(all_y, all_hat),
    }
    print(f"\n=== {title} ===")
    print(pd.DataFrame(fold_rows).to_string(index=False))
    print(
        f"Pooled AUROC={pooled['AUROC']:.4f} AUPRC={pooled['AUPRC']:.4f} MCC={pooled['MCC']:.4f}"
    )

    # also fit full for importances of new feats
    pipe.fit(X, y)
    imp = pd.Series(pipe.named_steps["classifier"].feature_importances_, index=feature_cols)
    new_imp = imp.reindex(NEW_FEATURES).dropna().sort_values(ascending=False)
    if len(new_imp):
        print("New-feature importances (full fit):")
        for k, v in new_imp.items():
            print(f"  {k}: {v:.4f}  (rank {int((imp > v).sum()) + 1}/{len(imp)})")
    pooled["importances_new"] = new_imp.to_dict()
    pooled["top10"] = imp.sort_values(ascending=False).head(10).to_dict()
    return pooled, fold_rows, imp


def main():
    print(f"loading {TRAIN}")
    df = pd.read_csv(TRAIN, sep="\t")
    optimal = pd.read_csv(OPT, sep="\t")["feature"].tolist()
    # drop header duplicate if present
    optimal = [c for c in optimal if c in df.columns]
    print(f"train rows={len(df)} optimal feats available={len(optimal)}")

    species_ids = set(df.loc[df["pr_cat"].isin(["TP", "FP"]), "species_cluster_id"].astype(int))
    annot_path = OUT_DIR / "liang_species_annot_features.tsv"
    if annot_path.exists():
        annot = pd.read_csv(annot_path, sep="\t")
        print(f"loaded cached annot features {annot_path}")
    else:
        annot = build_species_annot_features(species_ids)
        annot.to_csv(annot_path, sep="\t", index=False)
        print(f"wrote {annot_path}")

    df = df.merge(annot, on="species_cluster_id", how="left")
    for ft in NEW_FEATURES:
        df[ft] = df[ft].fillna(0.0)

    # prevalence in TP/FP
    sub = df[df["pr_cat"].isin(["TP", "FP"])]
    print("\nAnnot rates in TP vs FP (detection-level):")
    for ft in NEW_FEATURES:
        for cat in ["TP", "FP"]:
            s = sub[sub["pr_cat"] == cat][ft]
            print(f"  {ft} {cat}: {s.mean():.3f} (n={len(s)})")

    results = []
    configs = [
        ("baseline_86", optimal),
        ("plus_all_annot", optimal + NEW_FEATURES),
        ("plus_methylase", optimal + ["annot_has_methylase"]),
        ("plus_PHROG8", optimal + ["annot_has_PHROG8"]),
        ("plus_repr_exc", optimal + ["annot_has_repressor_excisionase"]),
        ("plus_methyl_PHROG8_repr", optimal + [
            "annot_has_methylase",
            "annot_has_PHROG8",
            "annot_has_repressor_excisionase",
        ]),
        ("plus_IS", optimal + ["annot_has_IS"]),
        ("annot_only", NEW_FEATURES),
    ]
    imp_store = {}
    for name, feats in configs:
        pooled, folds, imp = run_cv(df, feats, name)
        results.append(pooled)
        pd.DataFrame(folds).to_csv(OUT_DIR / f"folds_{name}.tsv", sep="\t", index=False)
        imp_store[name] = imp

    # univariate AUPRC of each new feature alone (as score: for inactive, invert if anti-correlated)
    print("\n=== Univariate signal of new features (AUPRC for inactive=TP) ===")
    y = (sub["pr_cat"] == "TP").astype(int)
    for ft in NEW_FEATURES:
        # score high = more inactive; if feature enriched in FP, use 1-x
        s = sub[ft].astype(float)
        # choose orientation by AUROC
        a1 = roc_auc_score(y, s)
        a0 = roc_auc_score(y, 1 - s)
        score = s if a1 >= a0 else 1 - s
        orient = "raw" if a1 >= a0 else "inverted"
        print(
            f"  {ft}: AUROC={max(a1,a0):.3f} AUPRC={average_precision_score(y, score):.3f} ({orient}) "
            f"P(TP|1)={y[s==1].mean():.3f} P(TP|0)={y[s==0].mean():.3f}"
        )

    summary = []
    base = next(r for r in results if r["title"] == "baseline_86")
    for r in results:
        summary.append(
            {
                "model": r["title"],
                "n_features": r["n_features"],
                "AUROC": r["AUROC"],
                "AUPRC": r["AUPRC"],
                "MCC": r["MCC"],
                "delta_AUROC": r["AUROC"] - base["AUROC"],
                "delta_AUPRC": r["AUPRC"] - base["AUPRC"],
                "delta_MCC": r["MCC"] - base["MCC"],
            }
        )
    summ_df = pd.DataFrame(summary)
    summ_df.to_csv(OUT_DIR / "annot_feature_ablation_summary.tsv", sep="\t", index=False)
    with open(OUT_DIR / "annot_feature_ablation_summary.json", "w") as f:
        json.dump({"results": results, "summary": summary}, f, indent=2, default=float)
    print("\n=== Summary vs baseline ===")
    print(summ_df.to_string(index=False))
    print(f"\nwrote {OUT_DIR}")


if __name__ == "__main__":
    main()
