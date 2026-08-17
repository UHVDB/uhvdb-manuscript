### Sylph taxonomy: ANI table + phage-to-host ratio (closest-to-1 primary)
# Primary VHR: closest-to-1 — among co-detected bacterial species in the predicted
# host genus (PHIST+CRISPR), pick the host with VHR nearest to 1.
# Cascade VHR is still computed as phage_host_ratio_cascade_df for comparison plots:
#   1) species_codetected — predicted host species co-detected in sample
#   2) genus_singleton — unique co-detected species within predicted host genus
#   3) family_singleton — unique co-detected species within predicted host family
# Host predictions: PHIST+CRISPR majority at species-cluster level (as in mmc_activity_analysis).

from pathlib import Path as _Path

FIG1_HOST = _Path('../figure_1/uhvdb_human_metag_results/uhvdb_2026-03-26-2')
PHIST_HOST = FIG1_HOST / 'uhvdb_phisthost.tsv.gz'
CRISPR_HOST = FIG1_HOST / 'uhvdb_crisprhost.tsv.gz'
# species_cluster_id must match sylph vSPECIES-* / figure_s18 metadata IDs
# (figure_1 uhvdb_species_info.cluster_id is a different ID space).
UHVDB_META = _Path('../figure_s18/uhvdb_v5_final_metadata_v2.tsv.gz')

uhvdb_genomovar_host = (
    pl.read_csv(PHIST_HOST, separator='\t')[
        ['uhvdb_id', 'total_connections', 'rank', 'agreement', 'consensus_taxonomy']
    ]
    .join(
        pl.read_csv(CRISPR_HOST, separator='\t')[
            ['uhvdb_id', 'total_connections', 'rank', 'agreement', 'top_taxonomy']
        ],
        on=['uhvdb_id', 'rank'],
        how='full',
        suffix='_crispr',
    )
    .fill_null(0)
    .with_columns(pl.col('consensus_taxonomy').str.replace(r'^(s|g|f)__', ''))
    .with_columns([
        pl.when(pl.col('consensus_taxonomy') == pl.col('top_taxonomy'))
        .then(pl.col('consensus_taxonomy'))
        .when(
            (pl.col('total_connections') * pl.col('agreement'))
            >= (pl.col('total_connections_crispr') * pl.col('agreement_crispr'))
        )
        .then(pl.col('consensus_taxonomy'))
        .otherwise(pl.col('top_taxonomy'))
        .alias('final_taxonomy'),
        pl.when(pl.col('rank').is_not_null()).then(pl.col('rank')).otherwise(pl.col('rank_crispr')).alias('final_rank'),
        pl.when(pl.col('uhvdb_id').is_not_null()).then(pl.col('uhvdb_id')).otherwise(pl.col('uhvdb_id_crispr')).alias('final_id'),
    ])
    .group_by('final_id')
    .agg([
        pl.col('final_taxonomy').filter(pl.col('final_rank') == 'species').first().alias('final_species'),
        pl.col('final_taxonomy').filter(pl.col('final_rank') == 'genus').first().alias('final_genus'),
        pl.col('final_taxonomy').filter(pl.col('final_rank') == 'family').first().alias('final_family'),
    ])
    .rename({'final_id': 'uhvdb_id'})
)

# Genomovar reps with correct species_cluster_id (matches sylph vSPECIES)
gv_in_species = (
    pl.read_csv(
        UHVDB_META, separator='\t',
        columns=['uhvdb_id', 'species_cluster_id', 'genomovar_rep', 'species_rep'],
    )
    .filter(pl.col('uhvdb_id') == pl.col('genomovar_rep'))
    .select(['uhvdb_id', 'species_cluster_id', 'species_rep'])
    .unique()
    .join(uhvdb_genomovar_host, on='uhvdb_id', how='left')
    .rename({'uhvdb_id': 'genomovar_rep'})
)


def majority_host(df: pl.DataFrame, col: str) -> pl.DataFrame:
    return (
        df
        .filter(pl.col(col).is_not_null())
        .group_by(['species_cluster_id', col])
        .agg(pl.len().alias('n'))
        .sort(['species_cluster_id', 'n', col], descending=[False, True, False])
        .unique('species_cluster_id', maintain_order=True)
        .select(['species_cluster_id', col])
    )


uhvdb_combinedhost = (
    gv_in_species
    .select(['species_cluster_id', 'species_rep'])
    .unique()
    .join(majority_host(gv_in_species, 'final_species'), on='species_cluster_id', how='left')
    .join(majority_host(gv_in_species, 'final_genus'), on='species_cluster_id', how='left')
    .join(majority_host(gv_in_species, 'final_family'), on='species_cluster_id', how='left')
    .select([
        'species_cluster_id',
        pl.col('species_rep').alias('uhvdb_id'),
        'final_species',
        'final_genus',
        'final_family',
    ])
)
print(
    'Host predictions: species',
    uhvdb_combinedhost.filter(pl.col('final_species').is_not_null()).height,
    '| genus',
    uhvdb_combinedhost.filter(pl.col('final_genus').is_not_null()).height,
    '| family',
    uhvdb_combinedhost.filter(pl.col('final_family').is_not_null()).height,
)

sylph_tax_results = []
virus_abund_lst = []
bac_abund_lst = []

for file in glob.glob('../figure_4/sylph_tax_results/*.sylphmpa'):
    sample_id = file.split('/')[-1].split('.syl')[0]
    df = pl.read_csv(
        file,
        separator='\t',
        skip_rows=1,
        new_columns=[
            'clade_name', 'taxonomic_abundance', 'relative_abundance',
            'ani', 'coverage', 'virus_host',
        ],
        null_values=['NA'],
        schema_overrides={
            'clade_name': pl.Utf8,
            'taxonomic_abundance': pl.Float64,
            'relative_abundance': pl.Float64,
            'ani': pl.Float64,
            'coverage': pl.Float64,
            'virus_host': pl.Utf8,
        },
    )

    virus_rows = (
        df
        .filter(pl.col('clade_name').str.starts_with('Viruses'))
        .filter(pl.col('clade_name').str.contains('t__'))
        .with_columns([
            pl.lit(sample_id).alias('sample_id'),
            pl.col('clade_name')
                .str.replace(r'.*vSPECIES-', '')
                .str.replace(r'\|.*', '')
                .cast(pl.Int64, strict=False)
                .alias('species_cluster_id'),
            pl.col('clade_name').str.extract(r't__(UHVDB-\d+)', 1).alias('uhvdb_id'),
        ])
        .rename({'taxonomic_abundance': 'virus_tax_abund'})
        .filter(pl.col('species_cluster_id').is_not_null())
        .filter(pl.col('virus_tax_abund') > 0)
    )
    sylph_tax_results.append(virus_rows.select(['sample_id', 'species_cluster_id', 'ani', 'uhvdb_id', 'virus_tax_abund']))
    virus_abund_lst.append(
        virus_rows
        .group_by(['sample_id', 'species_cluster_id'])
        .agg(pl.col('virus_tax_abund').max().alias('virus_tax_abund'))
    )

    bac_abund_lst.append(
        df
        .filter(pl.col('clade_name').str.starts_with('d__Bacteria'))
        .filter(pl.col('clade_name').str.contains('s__'))
        .with_columns([
            pl.lit(sample_id).alias('sample_id'),
            pl.col('clade_name').str.replace_all(';', '|'),
            pl.col('clade_name').str.extract(r's__([^|;]+)', 1).alias('species'),
            pl.col('clade_name').str.extract(r'g__([^|;]+)', 1).alias('genus'),
            pl.col('clade_name').str.extract(r'f__([^|;]+)', 1).alias('family'),
        ])
        .rename({'taxonomic_abundance': 'host_tax_abund'})
        .filter(pl.col('species').is_not_null())
        .filter(pl.col('host_tax_abund') > 0)
        .group_by(['sample_id', 'species', 'genus', 'family'])
        .agg(pl.col('host_tax_abund').max().alias('host_tax_abund'))
    )

sylph_tax_results_df = pl.concat(sylph_tax_results)
virus_abund = pl.concat(virus_abund_lst)
bac_sylph = pl.concat(bac_abund_lst)

viruses = (
    virus_abund
    .join(
        uhvdb_combinedhost.select([
            'species_cluster_id', 'final_species', 'final_genus', 'final_family',
        ]),
        on='species_cluster_id',
        how='left',
    )
    .with_columns([
        pl.col('final_genus').alias('host_genus'),
        pl.col('final_family').alias('host_family'),
    ])
)
print(f'Sylph virus rows (sample × species): {viruses.height:,}')
print(
    'with host species:', viruses.filter(pl.col('final_species').is_not_null()).height,
    '| genus:', viruses.filter(pl.col('host_genus').is_not_null()).height,
    '| family:', viruses.filter(pl.col('host_family').is_not_null()).height,
)

# 1) predicted host species co-detected
species_hits = (
    viruses
    .filter(pl.col('final_species').is_not_null() & (pl.col('virus_tax_abund') > 0))
    .join(
        bac_sylph.select(['sample_id', 'species', 'host_tax_abund']),
        left_on=['sample_id', 'final_species'],
        right_on=['sample_id', 'species'],
        how='inner',
    )
    .with_columns([
        pl.lit('species_codetected').alias('host_match_method'),
        pl.col('final_species').alias('matched_host_species'),
    ])
)
matched_keys = species_hits.select(['sample_id', 'species_cluster_id']).unique()

# 2) singleton species within predicted host genus
genus_singleton = (
    bac_sylph
    .filter(pl.col('genus').is_not_null())
    .group_by(['sample_id', 'genus'])
    .agg([
        pl.len().alias('n_species'),
        pl.col('species').first().alias('species'),
        pl.col('host_tax_abund').first().alias('host_tax_abund'),
    ])
    .filter(pl.col('n_species') == 1)
)
genus_hits = (
    viruses
    .join(matched_keys, on=['sample_id', 'species_cluster_id'], how='anti')
    .filter(pl.col('host_genus').is_not_null() & (pl.col('virus_tax_abund') > 0))
    .join(
        genus_singleton.select(['sample_id', 'genus', 'species', 'host_tax_abund']),
        left_on=['sample_id', 'host_genus'],
        right_on=['sample_id', 'genus'],
        how='inner',
    )
    .with_columns([
        pl.lit('genus_singleton').alias('host_match_method'),
        pl.col('species').alias('matched_host_species'),
    ])
)
matched_keys = pl.concat([
    matched_keys,
    genus_hits.select(['sample_id', 'species_cluster_id']),
]).unique()

# 3) singleton species within predicted host family
family_singleton = (
    bac_sylph
    .filter(pl.col('family').is_not_null())
    .group_by(['sample_id', 'family'])
    .agg([
        pl.len().alias('n_species'),
        pl.col('species').first().alias('species'),
        pl.col('host_tax_abund').first().alias('host_tax_abund'),
    ])
    .filter(pl.col('n_species') == 1)
)
family_hits = (
    viruses
    .join(matched_keys, on=['sample_id', 'species_cluster_id'], how='anti')
    .filter(pl.col('host_family').is_not_null() & (pl.col('virus_tax_abund') > 0))
    .join(
        family_singleton.select(['sample_id', 'family', 'species', 'host_tax_abund']),
        left_on=['sample_id', 'host_family'],
        right_on=['sample_id', 'family'],
        how='inner',
    )
    .with_columns([
        pl.lit('family_singleton').alias('host_match_method'),
        pl.col('species').alias('matched_host_species'),
    ])
)

vhr_cols = [
    'sample_id', 'species_cluster_id', 'virus_tax_abund', 'host_tax_abund',
    'final_species', 'host_genus', 'host_family',
    'matched_host_species', 'host_match_method',
]
# Cascade VHR kept for comparison with primary closest-to-1 method below.
phage_host_ratio_cascade_df = (
    pl.concat([
        species_hits.select(vhr_cols),
        genus_hits.select(vhr_cols),
        family_hits.select(vhr_cols),
    ])
    .with_columns((pl.col('virus_tax_abund') / pl.col('host_tax_abund')).alias('phage_host_ratio'))
    .filter(pl.col('host_tax_abund') > 0)
    .sort(['sample_id', 'species_cluster_id', 'virus_tax_abund'], descending=[False, False, True])
    .unique(subset=['sample_id', 'species_cluster_id'], keep='first')
)

print(f'Cascade VHR rows (matched host): {phage_host_ratio_cascade_df.height:,}')
print(phage_host_ratio_cascade_df.group_by('host_match_method').len().sort('host_match_method'))

# Primary VHR: closest-to-1 within predicted host genus (PHIST+CRISPR).
# Among all co-detected species in the predicted genus, keep the host whose
# virus/host abundance ratio is nearest to 1.
phage_host_ratio_df = (
    viruses
    .filter(pl.col('host_genus').is_not_null() & (pl.col('virus_tax_abund') > 0))
    .join(
        bac_sylph.select(['sample_id', 'species', 'genus', 'host_tax_abund']),
        left_on=['sample_id', 'host_genus'],
        right_on=['sample_id', 'genus'],
        how='inner',
    )
    .filter(pl.col('host_tax_abund') > 0)
    .with_columns(
        (pl.col('virus_tax_abund') / pl.col('host_tax_abund')).alias('phage_host_ratio')
    )
    .with_columns((pl.col('phage_host_ratio') - 1.0).abs().alias('dist_to_one'))
    .sort('dist_to_one')
    .unique(subset=['sample_id', 'species_cluster_id'], keep='first')
    .with_columns([
        pl.lit('closest_to_one').alias('host_match_method'),
        pl.col('species').alias('matched_host_species'),
    ])
    .select(vhr_cols + ['phage_host_ratio'])
)

print(f'Closest-to-1 VHR rows (matched host): {phage_host_ratio_df.height:,}')
print(phage_host_ratio_df.group_by('host_match_method').len().sort('host_match_method'))
