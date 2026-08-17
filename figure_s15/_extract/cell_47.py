# Classifier: inactive (TP) vs active (FP) among bulk-detected viruses.
# Feature notes:
#   - Includes bulk breadth / breadth_ratio / variance_ratio (bulk detection context).
#   - Excludes enriched breadth (label leakage) and host/ICTV ecology columns.
#   - med_prop_* gene-category proportions; SimpleImputer handles nulls (no gene-breadth).
#   - GroupKFold by sample group to reduce patient-level leakage.
COLS_TO_DROP = [
    # Identifiers & Targets
    "species_cluster_id",
    "sample_id",
    "group",
    "pr_cat",
    "is_active",
    "is_inactive",
    # Post-model analysis variable (must not be used to predict activity)
    "phage_host_ratio",
    # Leakage
    "breadth_enriched",
    "breadth_ratio_enriched",
    # Dropped to prevent ecological overfitting (Gut -> Airways)
    "most_common_ictv_class",
    "most_common_family_cluster_id",
    "most_common_host_taxonomy",
]


def run_activity_classifier(df_base, positive_pr_cat, model_title, target_col="target"):
    """Group K-fold RF classifier for bulk-detected viruses.

    positive_pr_cat: 'TP' = inactive (bulk, not enriched)
    """
    df = df_base[df_base["pr_cat"].isin(["TP", "FP"])].copy()
    df[target_col] = (df["pr_cat"] == positive_pr_cat).astype(int)

    n_pos = int(df[target_col].sum())
    n_neg = int((df[target_col] == 0).sum())
    print("=" * 60)
    print(model_title)
    print(f"Positive class ({positive_pr_cat}): {n_pos} | Negative class: {n_neg}")
    print("=" * 60)

    groups = df["group"]
    y = df[target_col]
    X = df.drop(columns=COLS_TO_DROP, errors="ignore")
    numeric_cols = X.columns.tolist()

    preprocessor = ColumnTransformer(transformers=[
        ("num", SimpleImputer(strategy="median"), numeric_cols),
    ])
    pipeline = Pipeline(steps=[
        ("preprocessor", preprocessor),
        ("classifier", RandomForestClassifier(
            n_estimators=100,
            class_weight="balanced",
            random_state=42,
            n_jobs=-1,
        )),
    ])

    n_splits = min(5, df["group"].nunique())
    gkf = GroupKFold(n_splits=n_splits)

    fold_results = []
    all_y_test, all_y_proba, all_y_pred = [], [], []
    oof_proba = pd.Series(np.nan, index=df.index, dtype=float)
    oof_pred = pd.Series(pd.NA, index=df.index, dtype="Int64")
    fold_roc_curves, fold_pr_curves = [], []

    print(f"Starting Group K-Fold Cross-Validation (k={n_splits})")

    for fold, (train_idx, test_idx) in enumerate(gkf.split(X, y, groups)):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        pipeline.fit(X_train, y_train)

        y_pred = pipeline.predict(X_test)
        y_proba = pipeline.predict_proba(X_test)[:, 1]
        test_index = df.index[test_idx]
        oof_proba.loc[test_index] = y_proba
        oof_pred.loc[test_index] = y_pred

        if len(np.unique(y_test)) > 1:
            auc = roc_auc_score(y_test, y_proba)
            auprc = average_precision_score(y_test, y_proba)
            mcc = matthews_corrcoef(y_test, y_pred)

            fpr, tpr, _ = roc_curve(y_test, y_proba)
            fold_roc_curves.append((fpr, tpr, auc))
            precision, recall, _ = precision_recall_curve(y_test, y_proba)
            fold_pr_curves.append((recall, precision, auprc))
        else:
            auc, auprc, mcc = np.nan, np.nan, np.nan

        fold_results.append({"fold": fold + 1, "AUROC": auc, "AUPRC": auprc, "MCC": mcc})
        all_y_test.extend(y_test.tolist())
        all_y_proba.extend(y_proba.tolist())
        all_y_pred.extend(y_pred.tolist())

    results_df = pd.DataFrame(fold_results).dropna()
    print(results_df.to_string(index=False))

    pooled_metrics = {}
    if len(np.unique(all_y_test)) > 1:
        pooled_metrics = {
            "AUROC": roc_auc_score(all_y_test, all_y_proba),
            "AUPRC": average_precision_score(all_y_test, all_y_proba),
            "MCC": matthews_corrcoef(all_y_test, all_y_pred),
        }
        print("\nPooled metrics across all folds:")
        print(f"  AUROC : {pooled_metrics['AUROC']:.3f}")
        print(f"  AUPRC : {pooled_metrics['AUPRC']:.3f}")
        print(f"  MCC   : {pooled_metrics['MCC']:.3f}")

    if len(fold_roc_curves) > 0:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        fig.suptitle(model_title, fontsize=14, fontweight="bold")

        ax = axes[0]
        for i, (fpr, tpr, auc) in enumerate(fold_roc_curves):
            ax.plot(fpr, tpr, lw=1.5, alpha=0.8, label=f"Fold {i + 1} (AUC = {auc:.3f})")

        pool_fpr, pool_tpr, _ = roc_curve(all_y_test, all_y_proba)
        pool_auc = roc_auc_score(all_y_test, all_y_proba)
        ax.plot(
            pool_fpr, pool_tpr, color="black", lw=2.5, linestyle="--",
            label=f"Pooled (AUC = {pool_auc:.3f})",
        )
        ax.plot([0, 1], [0, 1], color="grey", lw=1, linestyle=":")
        ax.set(xlabel="False Positive Rate", ylabel="True Positive Rate", title="ROC Curve (AUROC)")
        ax.legend(fontsize=9, loc="lower right")

        ax = axes[1]
        baseline = sum(all_y_test) / len(all_y_test)
        for i, (recall, precision, auprc) in enumerate(fold_pr_curves):
            ax.plot(recall, precision, lw=1.5, alpha=0.8, label=f"Fold {i + 1} (AP = {auprc:.3f})")

        pool_precision, pool_recall, _ = precision_recall_curve(all_y_test, all_y_proba)
        pool_ap = average_precision_score(all_y_test, all_y_proba)
        ax.plot(
            pool_recall, pool_precision, color="black", lw=2.5, linestyle="--",
            label=f"Pooled (AP = {pool_ap:.3f})",
        )
        ax.axhline(y=baseline, color="grey", lw=2, linestyle=":", label=f"Baseline ({baseline:.2f})")
        ax.axhline(y=0.95, color="black", lw=1.5, linestyle="--", label="Precision = 95%")
        ax.set(xlabel="Recall", ylabel="Precision", title="Precision-Recall Curve (AUPRC)")
        ax.set_ylim(0, 1)
        ax.legend(fontsize=9, loc="upper right")

        plt.tight_layout()
        plt.show()

    final_pipeline = Pipeline(steps=[
        ("preprocessor", preprocessor),
        ("classifier", RandomForestClassifier(
            n_estimators=100, class_weight="balanced", random_state=42, n_jobs=-1,
        )),
    ])
    final_pipeline.fit(X, y)

    X_transformed = final_pipeline.named_steps["preprocessor"].transform(X)
    X_transformed_df = pd.DataFrame(X_transformed, columns=numeric_cols)

    return {
        "df": df,
        "y": y,
        "X": X,
        "numeric_cols": numeric_cols,
        "results_df": results_df,
        "pooled_metrics": pooled_metrics,
        "fold_roc_curves": fold_roc_curves,
        "fold_pr_curves": fold_pr_curves,
        "all_y_test": all_y_test,
        "all_y_proba": all_y_proba,
        "oof_proba": oof_proba,
        "oof_pred": oof_pred,
        "final_pipeline": final_pipeline,
        "X_transformed_df": X_transformed_df,
        "positive_pr_cat": positive_pr_cat,
        "model_title": model_title,
    }


df_base = enrich_v_unenrich_final.to_pandas()

# Inactive viruses: detected in bulk but NOT enriched (pr_cat TP)
inactive_results = run_activity_classifier(
    df_base,
    positive_pr_cat="TP",
    model_title="Inactive Virus Model (positive = bulk, not enriched)",
    target_col="is_inactive",
)

# Names used by downstream SHAP cells
df = inactive_results["df"]
final_pipeline = inactive_results["final_pipeline"]
X_transformed_df = inactive_results["X_transformed_df"]
numeric_cols = inactive_results["numeric_cols"]
all_y_test = inactive_results["all_y_test"]
all_y_proba = inactive_results["all_y_proba"]
