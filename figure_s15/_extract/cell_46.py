# Feature table for classifier.
# med_prop_* = biological-group proportions (null if no gene-breadth data).
# Bulk breadth / breadth_ratio are intentional features for "inactive in bulk";
# enriched breadth columns are excluded as label leakage.
enrich_v_unenrich_final = (
    enrich_v_unenrich
    .with_columns([
        pl.when(pl.col('trimmed_mean') > 0)
            .then(pl.col('variance') / pl.col('trimmed_mean'))
            .otherwise(None)
            .alias('variance_ratio'),
    ])
    .select([
        'group', 'species_cluster_id', 'pr_cat',
        'breadth', 'breadth_ratio', 'variance_ratio', 'complete_count', 'high_quality_count',
        'med_virulent_score', 'med_n_hallmarks',
        'med_aai_id_af',
        'med_prop_capsid_packaging', 'med_prop_dna_metabolism',
        'med_prop_tail', 'med_prop_lysis', 'med_prop_connector',
        'med_prop_transcription', 'med_prop_integration',
        'med_prop_amg_host_takeover',
        'med_mcp_hallmark', 'med_portal_hallmark', 'med_terL_hallmark',
        # Retained for post-model analysis; explicitly excluded from classifier below.
        'phage_host_ratio',
        'ani',
    ])
    .sort(['group', 'species_cluster_id'])
)
enrich_v_unenrich_final.head(3)
