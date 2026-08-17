### Parse gene breadth results
# Biological-group *proportions* among genes with breadth > GENE_BREADTH_THRESHOLD.
# Requires shared definitions cell (with_bio_group_flags, BIO_GROUP_COLS).
# Hallmarks: presence if >= 1 matching gene among breadth-passing genes.

gene_breadth_lst = []

for file in glob.glob(
    '../figure_4/activity_profiling/liang_enriched_results/2026-04-03_outputs/'
    'referenceanalyze/coverm/*/bam/*.gene_coverage.tsv.gz'
):
    sample_id = file.split('/')[-1].split('.')[0]
    df = (
        pl.read_csv(file, separator='\t')
            .filter(pl.col('breadth') > GENE_BREADTH_THRESHOLD)
            .with_columns(pl.lit(sample_id).alias('sample_id'))
    )
    gene_breadth_lst.append(df)

gene_breadth_counts = (
    with_bio_group_flags(pl.concat(gene_breadth_lst))
        .group_by(['sample_id', 'genomovar_rep'])
        .agg([
            pl.len().alias('n_genes_covered'),
            pl.col('is_capsid_packaging').mean().alias('prop_capsid_packaging'),
            pl.col('is_dna_metabolism').mean().alias('prop_dna_metabolism'),
            pl.col('is_tail').mean().alias('prop_tail'),
            pl.col('is_lysis').mean().alias('prop_lysis'),
            pl.col('is_connector').mean().alias('prop_connector'),
            pl.col('is_transcription').mean().alias('prop_transcription'),
            pl.col('is_integration').mean().alias('prop_integration'),
            pl.col('is_amg_host_takeover').mean().alias('prop_amg_host_takeover'),
            pl.col('is_other').mean().alias('prop_other'),
            pl.col('is_unannotated').mean().alias('prop_unannotated'),
            (
                ((pl.col('pharokka_annot') == 'major head protein').sum() >= 1)
                | ((pl.col('phold_category') == 'major head protein').sum() >= 1)
                | ((pl.col('empathi_annot') == 'pvp|capsid|major_capsid').sum() >= 1)
            ).cast(pl.UInt32).alias('mcp_hallmark'),
            (
                ((pl.col('pharokka_annot') == 'terminase large subunit').sum() >= 1)
                | ((pl.col('phold_category') == 'terminase large subunit').sum() >= 1)
                | ((pl.col('empathi_annot') == 'DNA-associated|terminase|packaging_assembly').sum() >= 1)
            ).cast(pl.UInt32).alias('terl_hallmark'),
            (
                ((pl.col('pharokka_annot') == 'portal protein').sum() >= 1)
                | ((pl.col('phold_category') == 'portal protein').sum() >= 1)
                | ((pl.col('empathi_annot') == 'pvp|portal').sum() >= 1)
            ).cast(pl.UInt32).alias('portal_hallmark'),
        ])
        .with_columns(
            (
                pl.col('mcp_hallmark') + pl.col('terl_hallmark') + pl.col('portal_hallmark')
            ).alias('n_hallmarks')
        )
)
print(f'Gene breadth threshold: {GENE_BREADTH_THRESHOLD}')
print(f'Proportion columns: {[f"prop_{g}" for g in BIO_GROUP_COLS + ["unannotated"]]}')
gene_breadth_counts.select([
    'sample_id', 'genomovar_rep', 'n_genes_covered',
    'prop_capsid_packaging', 'prop_dna_metabolism', 'prop_tail', 'prop_lysis',
    'prop_connector', 'prop_transcription', 'prop_integration',
    'prop_amg_host_takeover', 'prop_other', 'prop_unannotated',
    'mcp_hallmark', 'terl_hallmark', 'portal_hallmark', 'n_hallmarks',
]).head(5)
