# TP = inactive (bulk only); FP = active (bulk + enriched); FN = enriched only
enrich_v_unenrich = (
    uhvdb_final_metadata_species
        .with_columns([
            pl.when((pl.col('breadth_ratio') > 0) & (pl.col('breadth_ratio_enriched') == 0))
                .then(pl.lit('TP'))
                .when((pl.col('breadth_ratio') == 0) & (pl.col('breadth_ratio_enriched') > 0))
                .then(pl.lit('FN'))
                .when((pl.col('breadth_ratio') > 0) & (pl.col('breadth_ratio_enriched') > 0))
                .then(pl.lit('FP'))
                .otherwise(pl.lit('TN'))
                .alias('pr_cat')
        ])
)


def _build_sylph_kmers_tables():
    """Load sylph .profile.tsv kmers_reassigned and aggregate to sample/group × species."""
    prof_glob = (
        '../figure_4/activity_profiling/liang_enriched_results/2026-04-03_outputs/'
        'referenceanalyze/sylph/*/*.profile.tsv'
    )
    lst = []
    for file in sorted(glob.glob(prof_glob)):
        sample_id = file.split('/')[-1].replace('.profile.tsv', '')
        df = (
            pl.read_csv(file, separator='\t')
            .filter(pl.col('Contig_name').str.starts_with('UHVDB-'))
            .with_columns([
                pl.lit(sample_id).alias('sample_id'),
                pl.col('Contig_name').alias('uhvdb_id'),
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
                'containment_den', 'Taxonomic_abundance',
            ])
        )
        lst.append(df)
    if not lst:
        raise FileNotFoundError(f'No sylph profile TSVs matched: {prof_glob}')
    id_map = (
        uhvdb_final_metadata
        .select(['uhvdb_id', 'species_cluster_id'])
        .unique('uhvdb_id')
        .with_columns(pl.col('species_cluster_id').cast(pl.Int64))
    )
    sp = (
        pl.concat(lst)
        .join(id_map, on='uhvdb_id', how='inner')
        .sort(['sample_id', 'species_cluster_id', 'Taxonomic_abundance'], descending=[False, False, True])
        .unique(['sample_id', 'species_cluster_id'], keep='first')
    )
    unenriched = (
        sp.filter(pl.col('sample_id').str.contains('_unenriched'))
        .select([
            'sample_id', 'species_cluster_id',
            'kmers_reassigned', 'kmers_reassigned_frac', 'containment_den',
        ])
    )
    enriched = (
        sp.filter(pl.col('sample_id').str.contains('_enriched'))
        .with_columns(pl.col('sample_id').str.replace('_enriched', '').alias('group'))
        .select([
            'group', 'species_cluster_id',
            pl.col('kmers_reassigned').alias('kmers_reassigned_enriched'),
            pl.col('kmers_reassigned_frac').alias('kmers_reassigned_frac_enriched'),
        ])
    )
    return unenriched, enriched


# Always (re)attach sylph kmers so this works even if the aggregate cell
# predated the kmers join or was not re-run.
_drop_kmers = [
    c for c in [
        'kmers_reassigned', 'kmers_reassigned_frac', 'containment_den',
        'kmers_reassigned_enriched', 'kmers_reassigned_frac_enriched',
    ]
    if c in enrich_v_unenrich.columns
]
if _drop_kmers:
    enrich_v_unenrich = enrich_v_unenrich.drop(_drop_kmers)

if 'sylph_kmers_by_species' in globals() and 'sylph_kmers_enriched_by_species' in globals():
    _km_un, _km_en = sylph_kmers_by_species, sylph_kmers_enriched_by_species
else:
    _km_un, _km_en = _build_sylph_kmers_tables()
    sylph_kmers_by_species, sylph_kmers_enriched_by_species = _km_un, _km_en

enrich_v_unenrich = (
    enrich_v_unenrich
    .with_columns(pl.col('species_cluster_id').cast(pl.Int64))
    .join(
        _km_un.with_columns(pl.col('species_cluster_id').cast(pl.Int64)),
        on=['sample_id', 'species_cluster_id'],
        how='left',
    )
    .join(
        _km_en.with_columns(pl.col('species_cluster_id').cast(pl.Int64)),
        on=['group', 'species_cluster_id'],
        how='left',
    )
)

print(enrich_v_unenrich['pr_cat'].value_counts().sort('pr_cat'))
print("Number of true positives:", enrich_v_unenrich.filter(pl.col('pr_cat') == 'TP').height)
print("Number of false negatives:", enrich_v_unenrich.filter(pl.col('pr_cat') == 'FN').height)
print("Number of false positives:", enrich_v_unenrich.filter(pl.col('pr_cat') == 'FP').height)
print(
    'With kmers_reassigned_frac:',
    enrich_v_unenrich.filter(pl.col('kmers_reassigned_frac').is_not_null()).height,
    '/', enrich_v_unenrich.height,
)
