# Feature table for classifier — optimal feature set
# (dump1 kitchen-sink + cm_dtr + n_struct_genes; best AUPRC from feature optimization).
# Bulk breadth / breadth_ratio are intentional; enriched breadth excluded (label leakage).
# Host / ICTV / source_db / VHR are excluded.

OPTIMAL_FEATURES = [
    'breadth',
    'breadth_ratio',
    'variance_ratio',
    'complete_count',
    'high_quality_count',
    'med_virulent_score',
    'med_n_hallmarks',
    'med_aai_id_af',
    'med_prop_capsid_packaging',
    'med_prop_dna_metabolism',
    'med_prop_tail',
    'med_prop_lysis',
    'med_prop_connector',
    'med_prop_transcription',
    'med_prop_integration',
    'med_prop_amg_host_takeover',
    'med_mcp_hallmark',
    'med_portal_hallmark',
    'med_terL_hallmark',
    'ani',
    'log_True_cov',
    'log_tax_abund',
    'log_contig_length',
    'log_n_genes',
    'completeness',
    'virus_score',
    'is_provirus',
    'host_genes',
    'viral_genes',
    'temperate',
    'annot_virulent',
    'phrog_integrases',
    'phrog_integration_excision',
    'empathi_integration',
    'is_integrated',
    'phrog_integrases_per_kb',
    'phrog_integration_excision_per_kb',
    'empathi_integration_per_kb',
    'num_lysis',
    'num_tail',
    'num_capsid',
    'num_proteins',
    'num_lysis_per_kb',
    'num_tail_per_kb',
    'num_capsid_per_kb',
    'num_proteins_per_kb',
    'annot_mcp_hallmark',
    'annot_terl_hallmark',
    'annot_portal_hallmark',
    'annot_n_hallmarks',
    'virus_hallmarks',
    'annot_n_hallmarks_per_kb',
    'virus_hallmarks_per_kb',
    'Adjusted_ANI',
    'Naive_ANI',
    'containment_frac',
    'Eff_lambda_ord',
    'ANI_interval_width',
    'Median_cov',
    'Mean_cov_geq1',
    'kmers_reassigned',
    'cov_evenness',
    'log_seq_abund',
    'gene_breadth_mean',
    'gene_breadth_std',
    'gene_breadth_cv',
    'frac_genes_breadth_ge_0_5',
    'frac_genes_breadth_ge_0_8',
    'gene_mean_depth_mean',
    'gene_mean_depth_std',
    'gene_occupancy',
    'annot_aai_id',
    'annot_aai_af',
    'aai_completeness',
    'aai_confidence_ord',
    'Score',
    'kmer_freq',
    'has_terminal_repeats',
    'is_DTR',
    'is_ITR',
    'is_topology_provirus',
    'proviral_frac',
    'True_cov_rank',
    'tax_abund_rank',
    'cm_dtr',
    'n_struct_genes'
]

_SYLPH_PROF_GLOB = (
    '../figure_4/activity_profiling/liang_enriched_results/2026-04-03_outputs/'
    'referenceanalyze/sylph/*/*.profile.tsv'
)
_GENE_COV_GLOB = (
    '../figure_4/activity_profiling/liang_enriched_results/2026-04-03_outputs/'
    'referenceanalyze/coverm/*/bam/*.gene_coverage.tsv.gz'
)


def _parse_containment_frac(s):
    if s is None or not isinstance(s, str) or '/' not in s:
        return None
    try:
        a, b = s.split('/', 1)
        a, b = float(a), float(b)
        return float(a / b) if b else None
    except Exception:
        return None


def _parse_ani_width(s: str):
    if s is None or not isinstance(s, str) or s.startswith('NA') or '-' not in s:
        return None
    try:
        lo, hi = s.split('-', 1)
        return float(hi) - float(lo)
    except Exception:
        return None


def _eff_lambda_ord(s: str):
    if s is None:
        return None
    v = {'LOW': 0.0, 'MEDIUM': 1.0, 'MED': 1.0, 'HIGH': 2.0}.get(str(s).upper())
    return v


# --- Species-rep metadata features (static genome annotations) ---
_meta_sp = (
    uhvdb_final_metadata
    .filter(pl.col('seq_name') == pl.col('seqhash_rep'))
    .unique('species_cluster_id', keep='first')
    .with_columns(pl.col('species_cluster_id').cast(pl.Int64))
    .with_columns([
        pl.col('virulent').alias('annot_virulent'),
        pl.col('temperate'),
        pl.col('phrog_integrases').cast(pl.Float64),
        pl.col('phrog_integration_excision').cast(pl.Float64),
        pl.col('empathi_integration').cast(pl.Float64),
        (pl.col('integration_status') == 'integrated').cast(pl.Float64).alias('is_integrated'),
        (pl.col('provirus').cast(pl.Utf8).str.to_lowercase().is_in(['yes', 'true', '1'])).cast(pl.Float64).alias('is_provirus'),
        pl.col('num_lysis').cast(pl.Float64), pl.col('num_tail').cast(pl.Float64), pl.col('num_capsid').cast(pl.Float64), pl.col('num_proteins').cast(pl.Float64),
        pl.col('mcp_hallmark').cast(pl.Float64).alias('annot_mcp_hallmark'),
        pl.col('terl_hallmark').cast(pl.Float64).alias('annot_terl_hallmark'),
        pl.col('portal_hallmark').cast(pl.Float64).alias('annot_portal_hallmark'),
        pl.col('n_hallmarks').cast(pl.Float64).alias('annot_n_hallmarks'),
        pl.col('virus_hallmarks').cast(pl.Float64),
        pl.col('aai_id').alias('annot_aai_id'),
        pl.col('aai_af').alias('annot_aai_af'),
        pl.col('aai_completeness'),
        pl.when(pl.col('aai_confidence') == 'low').then(0.0)
          .when(pl.col('aai_confidence') == 'medium').then(1.0)
          .when(pl.col('aai_confidence') == 'high').then(2.0)
          .otherwise(1.0)
          .alias('aai_confidence_ord'),
        pl.col('Score'),
        pl.col('kmer_freq'),
        pl.col('completeness'),
        pl.col('virus_score'),
        pl.col('host_genes'),
        pl.col('viral_genes'),
        pl.col('contig_length'),
        pl.col('n_genes'),
        pl.col('topology').is_in(['DTR', 'ITR']).cast(pl.Float64).alias('has_terminal_repeats'),
        (pl.col('topology') == 'DTR').cast(pl.Float64).alias('is_DTR'),
        (pl.col('topology') == 'ITR').cast(pl.Float64).alias('is_ITR'),
        (pl.col('topology') == 'Provirus').cast(pl.Float64).alias('is_topology_provirus'),
        (
            pl.col('proviral_length').cast(pl.Utf8).replace('NA', None).cast(pl.Float64, strict=False)
            / pl.col('contig_length')
        ).alias('proviral_frac'),
        pl.col('completeness_method').cast(pl.Utf8).str.contains('DTR').cast(pl.Float64).alias('cm_dtr'),
    ])
    .with_columns([
        (pl.col('num_lysis') / (pl.col('contig_length') / 1000)).alias('num_lysis_per_kb'),
        (pl.col('num_tail') / (pl.col('contig_length') / 1000)).alias('num_tail_per_kb'),
        (pl.col('num_capsid') / (pl.col('contig_length') / 1000)).alias('num_capsid_per_kb'),
        (pl.col('num_proteins') / (pl.col('contig_length') / 1000)).alias('num_proteins_per_kb'),
        (pl.col('phrog_integrases') / (pl.col('contig_length') / 1000)).alias('phrog_integrases_per_kb'),
        (pl.col('phrog_integration_excision') / (pl.col('contig_length') / 1000)).alias('phrog_integration_excision_per_kb'),
        (pl.col('empathi_integration') / (pl.col('contig_length') / 1000)).alias('empathi_integration_per_kb'),
        (pl.col('annot_n_hallmarks') / (pl.col('contig_length') / 1000)).alias('annot_n_hallmarks_per_kb'),
        (pl.col('virus_hallmarks') / (pl.col('contig_length') / 1000)).alias('virus_hallmarks_per_kb'),
        pl.col('contig_length').log1p().alias('log_contig_length'),
        pl.col('n_genes').log1p().alias('log_n_genes'),
    ])
    .select([
        'species_cluster_id',
        'annot_virulent', 'temperate', 'phrog_integrases', 'phrog_integration_excision', 'empathi_integration',
        'is_integrated', 'is_provirus',
        'num_lysis', 'num_tail', 'num_capsid', 'num_proteins',
        'num_lysis_per_kb', 'num_tail_per_kb', 'num_capsid_per_kb', 'num_proteins_per_kb',
        'phrog_integrases_per_kb', 'phrog_integration_excision_per_kb', 'empathi_integration_per_kb',
        'annot_mcp_hallmark', 'annot_terl_hallmark', 'annot_portal_hallmark', 'annot_n_hallmarks', 'virus_hallmarks',
        'annot_n_hallmarks_per_kb', 'virus_hallmarks_per_kb',
        'annot_aai_id', 'annot_aai_af', 'aai_completeness', 'aai_confidence_ord', 'Score', 'kmer_freq',
        'completeness', 'virus_score', 'host_genes', 'viral_genes',
        'log_contig_length', 'log_n_genes', 'n_genes',
        'has_terminal_repeats', 'is_DTR', 'is_ITR', 'is_topology_provirus', 'proviral_frac', 'cm_dtr',
    ])
)

# --- Extended Sylph profile features (unenriched / bulk training context) ---
_prof_lst = []
for _file in sorted(glob.glob(_SYLPH_PROF_GLOB)):
    _sid = Path(_file).name.replace('.profile.tsv', '')
    _df = (
        pl.read_csv(_file, separator='	')
        .filter(pl.col('Contig_name').str.starts_with('UHVDB-'))
        .with_columns([
            pl.lit(_sid).alias('sample_id'),
            pl.col('Contig_name').alias('uhvdb_id'),
            pl.col('True_cov').cast(pl.Float64),
            pl.col('Taxonomic_abundance').cast(pl.Float64).alias('tax_abund'),
            pl.col('Sequence_abundance').cast(pl.Float64).alias('seq_abund'),
            pl.col('Adjusted_ANI').cast(pl.Float64),
            pl.col('Naive_ANI').cast(pl.Float64),
            pl.col('Median_cov').cast(pl.Float64),
            pl.col('Mean_cov_geq1').cast(pl.Float64),
            pl.col('kmers_reassigned').cast(pl.Float64),
            pl.col('Containment_ind').cast(pl.Utf8),
            pl.col('Eff_lambda').cast(pl.Utf8),
            pl.col('ANI_5-95_percentile').cast(pl.Utf8).alias('ANI_5_95'),
        ])
    )
    _prof_lst.append(_df)

_id_map = (
    uhvdb_final_metadata
    .select(['uhvdb_id', 'species_cluster_id'])
    .unique('uhvdb_id')
    .with_columns(pl.col('species_cluster_id').cast(pl.Int64))
)

_prof_sp = (
    pl.concat(_prof_lst)
    .join(_id_map, on='uhvdb_id', how='inner')
    .sort(['sample_id', 'species_cluster_id', 'tax_abund'], descending=[False, False, True])
    .unique(['sample_id', 'species_cluster_id'], keep='first')
    .with_columns([
        pl.col('Containment_ind').map_elements(_parse_containment_frac, return_dtype=pl.Float64).alias('containment_frac'),
        (
            pl.when(pl.col('Eff_lambda').cast(pl.Utf8).str.to_uppercase() == 'LOW').then(0.0)
            .when(pl.col('Eff_lambda').cast(pl.Utf8).str.to_uppercase().is_in(['MEDIUM', 'MED'])).then(1.0)
            .when(pl.col('Eff_lambda').cast(pl.Utf8).str.to_uppercase() == 'HIGH').then(2.0)
            .otherwise(None)
            .alias('Eff_lambda_ord')
        ),
        pl.col('ANI_5_95').map_elements(_parse_ani_width, return_dtype=pl.Float64).alias('ANI_interval_width'),
        (pl.col('Median_cov') / pl.col('Mean_cov_geq1')).alias('cov_evenness'),
        pl.col('True_cov').log1p().alias('log_True_cov'),
        pl.col('tax_abund').fill_null(0).log1p().alias('log_tax_abund'),
        pl.col('seq_abund').fill_null(0).log1p().alias('log_seq_abund'),
    ])
)
_prof_sp = _prof_sp.with_columns([
    pl.col('True_cov').rank(method='average').over('sample_id').alias('_tc_rank'),
    pl.col('tax_abund').rank(method='average').over('sample_id').alias('_ta_rank'),
    pl.len().over('sample_id').alias('_n_in_sample'),
]).with_columns([
    (pl.col('_tc_rank') / pl.col('_n_in_sample')).alias('True_cov_rank'),
    (pl.col('_ta_rank') / pl.col('_n_in_sample')).alias('tax_abund_rank'),
]).select([
    'sample_id', 'species_cluster_id',
    'log_True_cov', 'log_tax_abund', 'log_seq_abund',
    'Adjusted_ANI', 'Naive_ANI', 'containment_frac', 'Eff_lambda_ord', 'ANI_interval_width',
    'Median_cov', 'Mean_cov_geq1', 'kmers_reassigned', 'cov_evenness',
    'True_cov_rank', 'tax_abund_rank',
])

# --- Gene coverage stats (all genes; includes n_struct_genes) ---
_struct = (
    (pl.col('pharokka_annot') == 'major head protein')
    | (pl.col('phold_category') == 'major head protein')
    | (pl.col('empathi_annot') == 'pvp|capsid|major_capsid')
    | (pl.col('pharokka_annot') == 'terminase large subunit')
    | (pl.col('phold_category') == 'terminase large subunit')
    | (pl.col('empathi_annot') == 'DNA-associated|terminase|packaging_assembly')
    | (pl.col('pharokka_annot') == 'portal protein')
    | (pl.col('phold_category') == 'portal protein')
    | (pl.col('empathi_annot') == 'pvp|portal')
)
_gene_lst = []
for _file in sorted(glob.glob(_GENE_COV_GLOB)):
    _sid = Path(_file).name.split('.')[0]
    _g = (
        pl.read_csv(_file, separator='	')
        .filter(pl.col('genomovar_rep').str.starts_with('UHVDB-'))
        .with_columns([
            pl.lit(_sid).alias('sample_id'),
            _struct.alias('is_struct'),
        ])
        .group_by(['sample_id', 'genomovar_rep'])
        .agg([
            pl.col('breadth').mean().alias('gene_breadth_mean'),
            pl.col('breadth').std().alias('gene_breadth_std'),
            (pl.col('breadth') >= 0.5).mean().alias('frac_genes_breadth_ge_0_5'),
            (pl.col('breadth') >= 0.8).mean().alias('frac_genes_breadth_ge_0_8'),
            pl.col('mean_depth').mean().alias('gene_mean_depth_mean'),
            pl.col('mean_depth').std().alias('gene_mean_depth_std'),
            pl.col('is_struct').sum().cast(pl.Float64).alias('n_struct_genes'),
            pl.len().alias('_n_genes_obs'),
            (pl.col('breadth') > GENE_BREADTH_THRESHOLD).sum().cast(pl.Float64).alias('n_genes_covered'),
        ])
        .with_columns(
            (pl.col('gene_breadth_std') / pl.col('gene_breadth_mean')).alias('gene_breadth_cv')
        )
    )
    _gene_lst.append(_g)

_gene_sp = (
    pl.concat(_gene_lst)
    .join(_id_map.rename({'uhvdb_id': 'genomovar_rep'}), on='genomovar_rep', how='inner')
    .sort(['sample_id', 'species_cluster_id', 'n_genes_covered'], descending=[False, False, True])
    .unique(['sample_id', 'species_cluster_id'], keep='first')
)

# Prefer species-rep n_genes from metadata for occupancy
_gene_sp = (
    _gene_sp
    .join(_meta_sp.select(['species_cluster_id', 'n_genes']), on='species_cluster_id', how='left')
    .with_columns(
        (pl.col('n_genes_covered') / pl.col('n_genes')).alias('gene_occupancy')
    )
    .select([
        'sample_id', 'species_cluster_id',
        'gene_breadth_mean', 'gene_breadth_std', 'gene_breadth_cv',
        'frac_genes_breadth_ge_0_5', 'frac_genes_breadth_ge_0_8',
        'gene_mean_depth_mean', 'gene_mean_depth_std', 'gene_occupancy', 'n_struct_genes',
    ])
)

# --- Assemble classifier feature table ---
enrich_v_unenrich_final = (
    enrich_v_unenrich
    .with_columns([
        pl.when(pl.col('trimmed_mean') > 0)
            .then(pl.col('variance') / pl.col('trimmed_mean'))
            .otherwise(None)
            .alias('variance_ratio'),
        pl.col('species_cluster_id').cast(pl.Int64),
    ])
    .join(_meta_sp, on='species_cluster_id', how='left')
    .join(_prof_sp, on=['sample_id', 'species_cluster_id'], how='left')
    .join(_gene_sp, on=['sample_id', 'species_cluster_id'], how='left')
)

_missing = [c for c in OPTIMAL_FEATURES if c not in enrich_v_unenrich_final.columns]
if _missing:
    raise KeyError(f'Optimal features missing from feature table: {_missing}')

enrich_v_unenrich_final = (
    enrich_v_unenrich_final
    .select(['group', 'species_cluster_id', 'pr_cat', 'phage_host_ratio'] + OPTIMAL_FEATURES)
    .sort(['group', 'species_cluster_id'])
)

print(f'Optimal feature set: {len(OPTIMAL_FEATURES)} features')
print(
    'Non-null coverage (bulk TP/FP) for key adds:',
    {
        c: enrich_v_unenrich_final.filter(pl.col('pr_cat').is_in(['TP', 'FP'])).select(c).drop_nulls().height
        for c in ['log_True_cov', 'containment_frac', 'gene_breadth_cv', 'cm_dtr', 'n_struct_genes', 'annot_aai_af']
    },
)
enrich_v_unenrich_final.head(3)
