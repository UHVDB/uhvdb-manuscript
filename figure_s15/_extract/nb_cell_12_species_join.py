### Characterize false negatives (enriched only; not in bulk)
# FN = breadth_ratio == 0 & breadth_ratio_enriched > 0
# Focus: how FN differ from FP (enriched detections that *are* also in bulk).

# Note: gene-breadth proportions are unavailable for FN rows (enriched-only detections);
# enriched CoverM outputs in this dataset lack *.gene_coverage.tsv.gz. Use ref_num_* instead.

fn = enrich_v_unenrich.filter(pl.col('pr_cat') == 'FN')
tp = enrich_v_unenrich.filter(pl.col('pr_cat') == 'TP')
fp = enrich_v_unenrich.filter(pl.col('pr_cat') == 'FP')

# If aggregate dropped enriched-only rows, rebuild and write back into enrich_v_unenrich
if fn.height == 0:
    fn = (
        uhvdb_final_metadata_species
        .filter((pl.col('breadth_ratio') == 0) & (pl.col('breadth_ratio_enriched') > 0))
        .with_columns(pl.lit('FN').alias('pr_cat'))
    )
    if fn.height > 0:
        enrich_v_unenrich = pl.concat(
            [enrich_v_unenrich.filter(pl.col('pr_cat') != 'FN'), fn],
            how='diagonal_relaxed',
        )
        print(f'Rebuilt FN and updated enrich_v_unenrich: {fn.height:,} rows')
    else:
        print(
            'WARNING: still 0 FN — re-run species-aggregate + pr_cat cells '
            '(check sample_id / contig_id coalesce for enriched-only rows).'
        )

fn = enrich_v_unenrich.filter(pl.col('pr_cat') == 'FN')
tp = enrich_v_unenrich.filter(pl.col('pr_cat') == 'TP')
fp = enrich_v_unenrich.filter(pl.col('pr_cat') == 'FP')

print(
    f"FN: {fn.height:,} rows | "
    f"{fn['species_cluster_id'].n_unique() if fn.height else 0:,} species | "
    f"{fn['group'].n_unique() if fn.height else 0:,} samples"
)
print(f"TP: {tp.height:,} rows | {tp['species_cluster_id'].n_unique():,} species")
print(f"FP: {fp.height:,} rows | {fp['species_cluster_id'].n_unique():,} species")

_meta_rep = (
    uhvdb_final_metadata
    .filter(pl.col('seq_name') == pl.col('seqhash_rep'))
    .select([
        'species_cluster_id', 'ictv_class', 'ictv_family', 'ictv_order',
        'host_lineage', 'family_cluster_id', 'contig_length', 'virulent', 'checkv_quality',
    ])
    .unique('species_cluster_id')
)


def _annotate(df):
    length_exprs = [pl.col('contig_length')]
    if 'genome_length' in df.columns:
        length_exprs = [pl.col('genome_length'), pl.col('contig_length')]
    return (
        df
        .join(_meta_rep, on='species_cluster_id', how='left')
        .with_columns([
            pl.col('host_lineage').str.replace_all(';', '|').str.extract(r'g__([^|;]+)', 1).alias('host_genus'),
            pl.col('host_lineage').str.replace_all(';', '|').str.extract(r'f__([^|;]+)', 1).alias('host_family'),
            pl.coalesce(length_exprs).alias('length_bp'),
        ])
    )


fn_annot = _annotate(fn)
tp_annot = _annotate(tp)
fp_annot = _annotate(fp)

if fn.height == 0:
    print('\nSkipping FN vs FP comparison (empty FN dataframe).')
else:
    # --- Numeric feature comparison: FN vs FP (TP shown for context) ---
    numeric_candidates = [
        'length_bp', 'breadth_enriched', 'breadth_ratio_enriched',
        'med_virulent_score', 'med_n_hallmarks', 'med_aai_id_af',
        'med_prop_capsid_packaging', 'med_prop_dna_metabolism',
        'med_prop_tail', 'med_prop_lysis', 'med_prop_connector',
        'med_prop_transcription', 'med_prop_integration',
        'med_mcp_hallmark', 'med_terL_hallmark', 'med_portal_hallmark',
        'phage_host_ratio', 'ani', 'complete_count', 'high_quality_count',
    ]
    numeric_cols = [
        c for c in numeric_candidates
        if c in fn_annot.columns and c in fp_annot.columns
    ]

    num_rows = []
    for col in numeric_cols:
        fn_vals = fn_annot[col].drop_nulls()
        fp_vals = fp_annot[col].drop_nulls()
        tp_vals = tp_annot[col].drop_nulls() if col in tp_annot.columns else None
        if fn_vals.len() < 3 or fp_vals.len() < 3:
            continue
        fn_mean = float(fn_vals.mean())
        fp_mean = float(fp_vals.mean())
        fn_std = float(fn_vals.std()) if fn_vals.len() > 1 else 0.0
        fp_std = float(fp_vals.std()) if fp_vals.len() > 1 else 0.0
        pooled = np.sqrt(
            ((fn_vals.len() - 1) * fn_std**2 + (fp_vals.len() - 1) * fp_std**2)
            / max(fn_vals.len() + fp_vals.len() - 2, 1)
        )
        d = (fn_mean - fp_mean) / pooled if pooled > 1e-12 else 0.0
        num_rows.append({
            'feature': col,
            'fn_mean': fn_mean,
            'fp_mean': fp_mean,
            'tp_mean': float(tp_vals.mean()) if tp_vals is not None and tp_vals.len() else None,
            'fn_median': float(fn_vals.median()),
            'fp_median': float(fp_vals.median()),
            'mean_diff_fn_minus_fp': fn_mean - fp_mean,
            'cohens_d_fn_vs_fp': d,
            'abs_cohens_d': abs(d),
        })

    fn_fp_num = pl.DataFrame(num_rows).sort('abs_cohens_d', descending=True)
    print('\n=== Numeric features: FN vs FP (ranked by |Cohen\'s d|) ===')
    print(fn_fp_num)

    # --- Presence / absence: % of viruses with feature == 0 (or null treated as absent for props) ---
    # This matches earlier FN findings (e.g. higher fraction with *no* integration genes).
    print(
        f"FN with gene-breadth props (med_prop_integration non-null): "
        f"{fn_annot.filter(pl.col('med_prop_integration').is_not_null()).height:,} / {fn_annot.height:,}"
    )
    if 'ref_num_integration' in fn_annot.columns:
        print(
            f"FN with ref_num_integration non-null: "
            f"{fn_annot.filter(pl.col('ref_num_integration').is_not_null()).height:,} / {fn_annot.height:,}"
        )

    presence_features = [
        # Gene-breadth proportions (0 = no genes in category among covered genes)
        'med_prop_integration', 'med_prop_lysis', 'med_prop_tail',
        'med_prop_capsid_packaging', 'med_prop_dna_metabolism',
        'med_prop_connector', 'med_prop_transcription',
        # Reference annotation count (earlier FN finding used this style of absence)
        'ref_num_integration',
        'med_n_hallmarks', 'med_mcp_hallmark', 'med_terL_hallmark', 'med_portal_hallmark',
    ]
    zero_rows = []
    for col in presence_features:
        if col not in fn_annot.columns or col not in fp_annot.columns:
            continue
        # Among rows with non-null gene-breadth features; null = no data (excluded)
        fn_obs = fn_annot.filter(pl.col(col).is_not_null())
        fp_obs = fp_annot.filter(pl.col(col).is_not_null())
        if fn_obs.height < 3 or fp_obs.height < 3:
            continue
        fn_zero = 100 * fn_obs.filter(pl.col(col) == 0).height / fn_obs.height
        fp_zero = 100 * fp_obs.filter(pl.col(col) == 0).height / fp_obs.height
        tp_obs = tp_annot.filter(pl.col(col).is_not_null()) if col in tp_annot.columns else None
        tp_zero = (
            100 * tp_obs.filter(pl.col(col) == 0).height / tp_obs.height
            if tp_obs is not None and tp_obs.height else None
        )
        zero_rows.append({
            'feature': col,
            'fn_pct_zero': fn_zero,
            'fp_pct_zero': fp_zero,
            'tp_pct_zero': tp_zero,
            'pct_point_diff_fn_minus_fp': fn_zero - fp_zero,
            'abs_pct_point_diff': abs(fn_zero - fp_zero),
            'fn_n': fn_obs.height,
            'fp_n': fp_obs.height,
        })

    fn_fp_zero = pl.DataFrame(zero_rows).sort('abs_pct_point_diff', descending=True)
    print('\n=== % viruses with feature == 0 (FN vs FP) ===')
    print('Positive pct_point_diff => more often *absent* in FN than FP')
    print(fn_fp_zero)

    MIN_ZERO_DIFF = 5.0
    top_zero = fn_fp_zero.filter(pl.col('abs_pct_point_diff') >= MIN_ZERO_DIFF).head(10)
    if top_zero.height == 0:
        top_zero = fn_fp_zero.head(min(6, fn_fp_zero.height))

    if top_zero.height > 0:
        plot_z = top_zero.sort('pct_point_diff_fn_minus_fp').to_pandas()
        fig, ax = plt.subplots(figsize=(8, max(3.5, 0.4 * len(plot_z))))
        colors = ['#c44e52' if v > 0 else '#4c72b0' for v in plot_z['pct_point_diff_fn_minus_fp']]
        y = range(len(plot_z))
        ax.barh(y, plot_z['pct_point_diff_fn_minus_fp'], color=colors)
        ax.set_yticks(list(y))
        ax.set_yticklabels([
            f"{r.feature}  (FN {r.fn_pct_zero:.0f}% vs FP {r.fp_pct_zero:.0f}% zero)"
            for r in plot_z.itertuples()
        ])
        ax.axvline(0, color='black', linewidth=0.8)
        ax.set_xlabel('Percentage-point difference in % zero (FN − FP)')
        ax.set_title('FN vs FP: gene/hallmark absence rates')
        plt.tight_layout()
        plt.show()


    TOP_N = 8
    MIN_D = 0.2
    top_num = fn_fp_num.filter(pl.col('abs_cohens_d') >= MIN_D).head(TOP_N)
    if top_num.height == 0:
        top_num = fn_fp_num.head(min(TOP_N, fn_fp_num.height))

    print(f'\nPlotting top {top_num.height} numeric features (|d| >= {MIN_D}, else top {TOP_N}):')
    print(top_num.select(['feature', 'fn_mean', 'fp_mean', 'mean_diff_fn_minus_fp', 'cohens_d_fn_vs_fp']))

    plot_df = top_num.sort('mean_diff_fn_minus_fp').to_pandas()
    fig, ax = plt.subplots(figsize=(8, max(3.5, 0.4 * len(plot_df))))
    colors = ['#c44e52' if v > 0 else '#4c72b0' for v in plot_df['mean_diff_fn_minus_fp']]
    y = range(len(plot_df))
    ax.barh(y, plot_df['mean_diff_fn_minus_fp'], color=colors)
    ax.set_yticks(list(y))
    ax.set_yticklabels([
        f"{r.feature}  (d={r.cohens_d_fn_vs_fp:.2f})"
        for r in plot_df.itertuples()
    ])
    ax.axvline(0, color='black', linewidth=0.8)
    ax.set_xlabel('Mean difference (FN − FP)')
    ax.set_title('FN vs FP: strongest numeric differences')
    plt.tight_layout()
    plt.show()

    n_dist = min(4, top_num.height)
    if n_dist > 0:
        fig, axes = plt.subplots(1, n_dist, figsize=(4 * n_dist, 4), squeeze=False)
        for ax, feat in zip(axes[0], top_num['feature'].to_list()[:n_dist]):
            parts = []
            for label, df in [('FN', fn_annot), ('FP', fp_annot), ('TP', tp_annot)]:
                vals = df[feat].drop_nulls().to_list() if feat in df.columns else []
                parts.append(vals)
            ax.boxplot(parts, tick_labels=['FN', 'FP', 'TP'], showfliers=False)
            ax.set_title(feat.replace('med_', ''), fontsize=10)
            ax.set_ylabel('Value')
        fig.suptitle('FN vs FP vs TP distributions (top features)', y=1.02)
        plt.tight_layout()
        plt.show()

    # --- Categorical composition: FN vs FP (% point differences) ---
    cat_cols = [
        ('host_genus', 'Host genus'),
        ('host_family', 'Host family'),
        ('ictv_class', 'ICTV class'),
        ('ictv_family', 'ICTV family'),
        ('checkv_quality', 'CheckV quality'),
    ]
    cat_diff_rows = []
    for col, title in cat_cols:
        if col not in fn_annot.columns:
            continue
        fn_pct = (
            fn_annot.group_by(col).agg(pl.len().alias('n'))
            .with_columns((100 * pl.col('n') / fn_annot.height).alias('pct_fn'))
            .select([col, 'pct_fn'])
        )
        fp_pct = (
            fp_annot.group_by(col).agg(pl.len().alias('n'))
            .with_columns((100 * pl.col('n') / max(fp_annot.height, 1)).alias('pct_fp'))
            .select([col, 'pct_fp'])
        )
        merged = (
            fn_pct.join(fp_pct, on=col, how='full', coalesce=True)
            .with_columns([
                pl.col('pct_fn').fill_null(0.0),
                pl.col('pct_fp').fill_null(0.0),
            ])
            .with_columns(
                (pl.col('pct_fn') - pl.col('pct_fp')).alias('pct_diff_fn_minus_fp'),
                (pl.col('pct_fn') - pl.col('pct_fp')).abs().alias('abs_pct_diff'),
                pl.lit(title).alias('category_type'),
                pl.col(col).cast(pl.Utf8).alias('level'),
            )
            .sort('abs_pct_diff', descending=True)
        )
        cat_diff_rows.append(merged.select([
            'category_type', 'level', 'pct_fn', 'pct_fp', 'pct_diff_fn_minus_fp', 'abs_pct_diff',
        ]))

    if cat_diff_rows:
        cat_diff = pl.concat(cat_diff_rows).sort('abs_pct_diff', descending=True)
        MIN_PCT = 5.0
        top_cat = cat_diff.filter(pl.col('abs_pct_diff') >= MIN_PCT).head(12)
        if top_cat.height == 0:
            top_cat = cat_diff.head(8)

        print('\n=== Categorical composition: FN vs FP (largest % point diffs) ===')
        print(top_cat)

        plot_cat = top_cat.sort('pct_diff_fn_minus_fp').to_pandas()
        fig, ax = plt.subplots(figsize=(9, max(3.5, 0.35 * len(plot_cat))))
        colors = ['#c44e52' if v > 0 else '#4c72b0' for v in plot_cat['pct_diff_fn_minus_fp']]
        y = range(len(plot_cat))
        ax.barh(y, plot_cat['pct_diff_fn_minus_fp'], color=colors)
        ax.set_yticks(list(y))
        ax.set_yticklabels([
            f"{r.category_type}: {r.level}" for r in plot_cat.itertuples()
        ])
        ax.axvline(0, color='black', linewidth=0.8)
        ax.set_xlabel('Percentage-point difference (FN − FP)')
        ax.set_title('FN vs FP: strongest categorical enrichment differences')
        plt.tight_layout()
        plt.show()

    print('\n=== Length bins (unique species): FN vs FP ===')

    def _len_bins(df):
        return (
            df.unique('species_cluster_id')
            .with_columns(
                pl.when(pl.col('length_bp') < 15_000).then(pl.lit('<15kb'))
                .when(pl.col('length_bp') < 30_000).then(pl.lit('15-30kb'))
                .when(pl.col('length_bp') < 50_000).then(pl.lit('30-50kb'))
                .when(pl.col('length_bp') < 100_000).then(pl.lit('50-100kb'))
                .otherwise(pl.lit('≥100kb'))
                .alias('len_bin')
            )
            .group_by('len_bin')
            .agg(pl.len().alias('n'))
            .with_columns((100 * pl.col('n') / pl.col('n').sum()).alias('pct'))
        )

    len_cmp = (
        _len_bins(fn_annot).rename({'n': 'n_fn', 'pct': 'pct_fn'})
        .join(
            _len_bins(fp_annot).rename({'n': 'n_fp', 'pct': 'pct_fp'}),
            on='len_bin', how='full', coalesce=True,
        )
        .with_columns([
            pl.col('pct_fn').fill_null(0.0),
            pl.col('pct_fp').fill_null(0.0),
        ])
        .with_columns(
            (pl.col('pct_fn') - pl.col('pct_fp')).alias('pct_diff')
        )
        .sort('pct_diff', descending=True)
    )
    print(len_cmp)
