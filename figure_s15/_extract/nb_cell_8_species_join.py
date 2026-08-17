# aggregate to species level (species representatives only; no genomovar median)
# Gene-category features: med_prop_* from gene_breadth_counts (unenriched samples only —
# enriched CoverM dirs in this dataset have no *.gene_coverage.tsv.gz outputs).
# Hallmarks / med_prop_* stay null when gene-breadth is missing (not filled as 0).

# Sylph genome-level profiles (kmers_reassigned) for unenriched + enriched libraries.
# % kmers reassigned = kmers_reassigned / containment denominator (Containment_ind "a/b").
_SYLPH_PROF_GLOB = (
    '../figure_4/activity_profiling/liang_enriched_results/2026-04-03_outputs/'
    'referenceanalyze/sylph/*/*.profile.tsv'
)
_sylph_prof_lst = []
for _file in sorted(glob.glob(_SYLPH_PROF_GLOB)):
    _sample_id = _file.split('/')[-1].replace('.profile.tsv', '')
    _df = (
        pl.read_csv(_file, separator='\t')
        .filter(pl.col('Contig_name').str.starts_with('UHVDB-'))
        .with_columns([
            pl.lit(_sample_id).alias('sample_id'),
            pl.col('Contig_name').alias('uhvdb_id'),
            pl.col('Containment_ind').str.split('/').list.get(0).cast(pl.Float64, strict=False).alias('containment_num'),
            pl.col('Containment_ind').str.split('/').list.get(1).cast(pl.Float64, strict=False).alias('containment_den'),
            pl.col('kmers_reassigned').cast(pl.Float64, strict=False),
        ])
        .with_columns(
            pl.when(pl.col('containment_den') > 0)
            .then(pl.col('kmers_reassigned') / pl.col('containment_den'))
            .otherwise(None)
            .alias('kmers_reassigned_frac')
        )
        .select([
            'sample_id', 'uhvdb_id', 'kmers_reassigned', 'kmers_reassigned_frac',
            'containment_den', 'Taxonomic_abundance', 'Adjusted_ANI',
        ])
    )
    _sylph_prof_lst.append(_df)

sylph_profile_df = pl.concat(_sylph_prof_lst) if _sylph_prof_lst else pl.DataFrame()
print(f'Sylph profile virus rows: {sylph_profile_df.height:,} across {sylph_profile_df["sample_id"].n_unique() if sylph_profile_df.height else 0:,} samples')

# Map UHVDB contig -> species_cluster_id (profiles are species-rep oriented)
_sylph_id_map = (
    uhvdb_final_metadata
    .select(['uhvdb_id', 'species_cluster_id'])
    .unique('uhvdb_id')
    .with_columns(pl.col('species_cluster_id').cast(pl.Int64))
)
_sylph_sp = (
    sylph_profile_df
    .join(_sylph_id_map, on='uhvdb_id', how='inner')
    .sort(['sample_id', 'species_cluster_id', 'Taxonomic_abundance'], descending=[False, False, True])
    .unique(['sample_id', 'species_cluster_id'], keep='first')
)

# Primary: unenriched profiles (match enrich_v_unenrich sample_id for TP/FP)
sylph_kmers_by_species = (
    _sylph_sp
    .filter(pl.col('sample_id').str.contains('_unenriched'))
    .select([
        'sample_id', 'species_cluster_id',
        'kmers_reassigned', 'kmers_reassigned_frac', 'containment_den',
    ])
)
# Enriched profiles keyed by group for optional join
sylph_kmers_enriched_by_species = (
    _sylph_sp
    .filter(pl.col('sample_id').str.contains('_enriched'))
    .with_columns(pl.col('sample_id').str.replace('_enriched', '').alias('group'))
    .select([
        'group', 'species_cluster_id',
        pl.col('kmers_reassigned').alias('kmers_reassigned_enriched'),
        pl.col('kmers_reassigned_frac').alias('kmers_reassigned_frac_enriched'),
    ])
)
print(f'Sylph kmers unenriched sample×species: {sylph_kmers_by_species.height:,}')
print(f'Sylph kmers enriched group×species: {sylph_kmers_enriched_by_species.height:,}')

_joined = (
    coverm_df
        .filter(pl.col('sample_id').str.contains('_unenriched'))
        .join(
            coverm_df.filter(pl.col('sample_id').str.contains('_enriched')),
            on=['group', 'contig_id'], suffix='_enriched', how='full'
        )
)

_coalesce_exprs = []
if 'contig_id_enriched' in _joined.columns:
    _coalesce_exprs.append(
        pl.coalesce([pl.col('contig_id'), pl.col('contig_id_enriched')]).alias('contig_id')
    )
if 'group_enriched' in _joined.columns:
    _coalesce_exprs.append(
        pl.coalesce([pl.col('group'), pl.col('group_enriched')]).alias('group')
    )

_gb = gene_breadth_counts.rename({'genomovar_rep': 'contig_id'})

uhvdb_final_metadata_species = (
    (_joined.with_columns(_coalesce_exprs) if _coalesce_exprs else _joined)
        .filter(pl.col('contig_id').is_not_null())
        .join(uhvdb_final_metadata, left_on='contig_id', right_on='uhvdb_id', how='left')
        .filter(pl.col('seq_name') == pl.col('seqhash_rep'))
        .drop([
            'num_capsid', 'num_tail', 'num_lysis',
            'mcp_hallmark', 'terl_hallmark', 'portal_hallmark',
            'n_hallmarks',
        ], strict=False)
        # Bulk (unenriched) gene-breadth only — enriched gene coverage was not generated
        .join(_gb, on=['contig_id', 'sample_id'], how='left')
        .with_columns([
            # Analysis sample_id: bulk if present, else enriched (FN)
            pl.coalesce([pl.col('sample_id'), pl.col('sample_id_enriched')]).alias('sample_id'),
            # Reference integration signal (metadata annotations; num_integration does not exist)
            (
                pl.col('phrog_integrases').fill_null(0)
                + pl.col('phrog_integration_excision').fill_null(0)
                + pl.col('empathi_integration').fill_null(0)
            ).alias('ref_num_integration'),
        ])
        .with_columns([
            (pl.col('checkv_quality') == 'Complete').cast(pl.Int64).alias('complete_count'),
            (pl.col('checkv_quality') == 'High-quality').cast(pl.Int64).alias('high_quality_count'),
            pl.col('n_hallmarks').alias('med_n_hallmarks'),
            ((pl.col('aai_id') / 100) * pl.col('aai_af')).alias('med_aai_id_af'),
            pl.col('ictv_class').alias('most_common_ictv_class'),
            pl.col('family_cluster_id').alias('most_common_family_cluster_id'),
            pl.col('host_lineage').alias('most_common_host_taxonomy'),
            pl.col('phist_genus_connections').alias('med_phist_genus_connections'),
            pl.col('phist_species_connections').alias('med_phist_species_connections'),
            pl.col('phist_family_connections').alias('med_phist_family_connections'),
            pl.col('crispr_genus_connections').alias('med_crispr_genus_connections'),
            pl.col('crispr_species_connections').alias('med_crispr_species_connections'),
            pl.col('crispr_family_connections').alias('med_crispr_family_connections'),
            pl.col('virulent').alias('med_virulent_score'),
            pl.col('prop_capsid_packaging').alias('med_prop_capsid_packaging'),
            pl.col('prop_dna_metabolism').alias('med_prop_dna_metabolism'),
            pl.col('prop_tail').alias('med_prop_tail'),
            pl.col('prop_lysis').alias('med_prop_lysis'),
            pl.col('prop_connector').alias('med_prop_connector'),
            pl.col('prop_transcription').alias('med_prop_transcription'),
            pl.col('prop_integration').alias('med_prop_integration'),
            pl.col('prop_amg_host_takeover').alias('med_prop_amg_host_takeover'),
            pl.col('prop_other').alias('med_prop_other'),
            pl.col('prop_unannotated').alias('med_prop_unannotated'),
            (pl.col('num_uniprot_ips') / pl.col('num_proteins')).alias('mean_proportion_uniprot_ips'),
            pl.col('mcp_hallmark').alias('med_mcp_hallmark'),
            pl.col('terl_hallmark').alias('med_terL_hallmark'),
            pl.col('portal_hallmark').alias('med_portal_hallmark'),
            pl.col('length').alias('genome_length'),
        ])
        .select([
            'species_cluster_id', 'sample_id', 'group',
            'complete_count', 'high_quality_count',
            'med_n_hallmarks', 'med_aai_id_af',
            'most_common_ictv_class', 'most_common_family_cluster_id', 'most_common_host_taxonomy',
            'med_phist_genus_connections', 'med_phist_species_connections', 'med_phist_family_connections',
            'med_crispr_genus_connections', 'med_crispr_species_connections', 'med_crispr_family_connections',
            'med_virulent_score', 'mean_proportion_uniprot_ips',
            'med_prop_capsid_packaging', 'med_prop_dna_metabolism',
            'med_prop_tail', 'med_prop_lysis', 'med_prop_connector',
            'med_prop_transcription', 'med_prop_integration',
            'med_prop_amg_host_takeover', 'med_prop_other', 'med_prop_unannotated',
            'ref_num_integration',
            'med_mcp_hallmark', 'med_terL_hallmark', 'med_portal_hallmark',
            'breadth', 'breadth_enriched', 'breadth_ratio', 'breadth_ratio_enriched',
            'variance', 'trimmed_mean', 'genome_length',
        ])
        .with_columns(pl.col('species_cluster_id').cast(pl.Int64))
        .join(
            phage_host_ratio_df
                .select(['sample_id', 'species_cluster_id', 'phage_host_ratio', 'host_match_method', 'matched_host_species'])
                .with_columns(pl.col('species_cluster_id').cast(pl.Int64)),
            on=['sample_id', 'species_cluster_id'],
            how='left',
        )
        .join(
            sylph_tax_results_df
                .select(['sample_id', 'species_cluster_id', 'ani'])
                .unique(['sample_id', 'species_cluster_id'])
                .with_columns(pl.col('species_cluster_id').cast(pl.Int64)),
            on=['sample_id', 'species_cluster_id'],
            how='left',
        )
        .join(
            sylph_kmers_by_species
                .with_columns(pl.col('species_cluster_id').cast(pl.Int64)),
            on=['sample_id', 'species_cluster_id'],
            how='left',
        )
        .join(
            sylph_kmers_enriched_by_species
                .with_columns(pl.col('species_cluster_id').cast(pl.Int64)),
            on=['group', 'species_cluster_id'],
            how='left',
        )
        .unique(['species_cluster_id', 'sample_id'])
        # Fill only coverm coverage/nulls from the full join — never gene-breadth props/hallmarks.
        .with_columns([
            pl.col(c).fill_null(0.0)
            for c in [
                'breadth', 'breadth_enriched', 'breadth_ratio', 'breadth_ratio_enriched',
                'variance', 'trimmed_mean',
                'complete_count', 'high_quality_count',
            ]
        ])
        .with_columns(pl.col('species_cluster_id').cast(pl.Int64))
)

_fn = uhvdb_final_metadata_species.filter(
    (pl.col('breadth_ratio') == 0) & (pl.col('breadth_ratio_enriched') > 0)
)
print('FN candidates (enriched only):', _fn.height)
print(
    'FN with non-null med_prop_integration (expected 0; no enriched gene coverage):',
    _fn.filter(pl.col('med_prop_integration').is_not_null()).height,
)
print(
    'FN with non-null ref_num_integration:',
    _fn.filter(pl.col('ref_num_integration').is_not_null()).height,
)
uhvdb_final_metadata_species.head(1)
