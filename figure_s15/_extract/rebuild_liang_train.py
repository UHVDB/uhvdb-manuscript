#!/usr/bin/env python3
"""Rebuild Liang enrich_v_unenrich_final and test annotation feature additions."""
import os
from pathlib import Path
os.chdir(Path(__file__).resolve().parents[1])
print("cwd", os.getcwd())

# ===== notebook cell 1 =====
### load packages
import math
import glob
from pathlib import Path
import os
import random

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
import seaborn as sns
import shap
from scipy import stats as scipy_stats
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import (
    average_precision_score,
    matthews_corrcoef,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline

plt.rcParams.update({'font.size': 14})

# Set all possible seeds
np.random.seed(1)
random.seed(1)
os.environ['PYTHONHASHSEED'] = '1'


# ===== notebook cell 3 =====
### load coverm results to identify viruses detected in bulk and enriched
df_lst = []

for file in glob.glob('../figure_4/activity_profiling/liang_enriched_results/2026-04-03_outputs/referenceanalyze/coverm/*/*.coverm.tsv.gz'):
    sample_id = file.split('/')[-1].split('.')[0]
    group = sample_id.rsplit('_', 1)[0]
    df = (
        pl.read_csv(file, separator='\t', new_columns=['contig_id', 'trimmed_mean', 'mean', 'variance', 'covered_bases', 'length'])
            .with_columns([
                pl.lit(sample_id).alias('sample_id'),
                pl.lit(group).alias('group'),
            ])
            .with_columns([
                (pl.col('covered_bases') / pl.col('length')).alias('breadth'),
                (1 - math.e**(-0.833 * pl.col('mean'))).alias('expected_breadth'),
            ])
            .with_columns([
                # Avoid divide-by-near-zero when expected breadth is tiny
                pl.when(pl.col('expected_breadth') > 1e-6)
                    .then(pl.col('breadth') / pl.col('expected_breadth'))
                    .otherwise(None)
                    .alias('breadth_ratio'),
            ])
    )
    df_lst.append(df)

coverm_df = pl.concat(df_lst)


# ===== notebook cell 4 =====
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


# ===== notebook cell 5 =====
### Gene-breadth biological groups (shared definitions)
# Proportion metric:
#   prop_group = (# genes with breadth > threshold in group)
#              / (# genes with breadth > threshold)
# Membership (tightened for manuscript):
#   Pharokka: exact category strings
#   phold: category/annot substring matches (phold is free-text)
#   Empathi: exact token0 / token1 only (no broad annot.contains)
# Run this cell before parse / discrimination.

GENE_BREADTH_THRESHOLD = 0.8

BIO_GROUP_COLS = [
    'capsid_packaging',
    'dna_metabolism',
    'tail',
    'lysis',
    'connector',
    'transcription',
    'integration',
    'amg_host_takeover',
    'other',
]

# Values that do not provide a functional annotation. This is intentionally
# strict: a gene is unannotated only when Pharokka, Phold, and Empathi are all
# uninformative, not merely when it falls outside BIO_GROUP_COLS.
UNINFORMATIVE_ANNOTATION_RE = (
    r'(?i)^\s*(?:hypothetical(?: protein)?|uncharacteri[sz]ed(?: protein)?|'
    r'unknown(?: protein| function)?|no annotation|no[_ -]?phrog|none|na|n/a|-)?\s*$'
)


def with_bio_group_flags(df: pl.DataFrame) -> pl.DataFrame:
    """Add biological-group flags and a strict all-sources-unannotated flag."""
    cols = df.columns
    prep = []
    if 'pharokka_category' in cols:
        prep.append(pl.col('pharokka_category').fill_null(''))
    else:
        prep.append(pl.lit('').alias('pharokka_category'))
    if 'phold_category' in cols:
        prep.append(pl.col('phold_category').fill_null(''))
    else:
        prep.append(pl.lit('').alias('phold_category'))
    if 'empathi_annot' in cols:
        prep.append(pl.col('empathi_annot').fill_null(''))
    else:
        prep.append(pl.lit('').alias('empathi_annot'))
    if 'pharokka_annot' in cols:
        prep.append(pl.col('pharokka_annot').fill_null(''))
    else:
        prep.append(pl.lit('').alias('pharokka_annot'))

    out = df.with_columns(prep).with_columns([
        pl.col('empathi_annot')
            .str.split('|')
            .list.get(0, null_on_oob=True)
            .fill_null('')
            .alias('empathi_token0'),
        pl.col('empathi_annot')
            .str.split('|')
            .list.get(1, null_on_oob=True)
            .fill_null('')
            .alias('empathi_token1'),
    ])
    return out.with_columns([
        (
            (pl.col('pharokka_category') == 'head and packaging')
            | pl.col('phold_category').str.contains('(?i)head|capsid|portal|terminase')
            | (pl.col('empathi_token0') == 'pvp')
            | (pl.col('empathi_token0') == 'packaging_assembly')
            | pl.col('empathi_token1').is_in([
                'capsid', 'terminase', 'portal', 'head-tail_joining',
            ])
        ).alias('is_capsid_packaging'),
        (
            (pl.col('pharokka_category') == 'DNA, RNA and nucleotide metabolism')
            | (pl.col('empathi_token0') == 'DNA-associated')
            | (pl.col('empathi_token0') == 'RNA-associated')
            | pl.col('empathi_token1').is_in([
                'nuclease', 'annealing', 'DNA_polymerase', 'helicase',
            ])
        ).alias('is_dna_metabolism'),
        (
            (pl.col('pharokka_category') == 'tail')
            | pl.col('phold_category').str.contains('(?i)\\btail\\b')
            | (pl.col('empathi_token1') == 'tail')
        ).alias('is_tail'),
        (
            (pl.col('pharokka_category') == 'lysis')
            | pl.col('phold_category').str.contains('(?i)lysis|holin|endolysin')
            | (pl.col('empathi_token0') == 'lysis')
            | (pl.col('empathi_token0') == 'cell_wall_depolymerase')
            | pl.col('empathi_token1').is_in(['lysis', 'holin'])
        ).alias('is_lysis'),
        (
            (pl.col('pharokka_category') == 'connector')
            | pl.col('phold_category').str.contains('(?i)connector|head-tail|head–tail')
        ).alias('is_connector'),
        (
            (pl.col('pharokka_category') == 'transcription regulation')
            | (pl.col('empathi_token0') == 'transcriptional_regulator')
            | (pl.col('empathi_token1') == 'transcriptional_regulator')
        ).alias('is_transcription'),
        (
            (pl.col('pharokka_category') == 'integration and excision')
            | pl.col('phold_category').str.contains('(?i)integrase|integration')
            | (pl.col('empathi_token1') == 'integration')
        ).alias('is_integration'),
        (
            (pl.col('pharokka_category') == 'moron, auxiliary metabolic gene and host takeover')
            | pl.col('phold_category').str.contains(
                '(?i)auxiliary metabolic|host takeover|moron|anti[- ]?defen[cs]e'
            )
            | pl.col('empathi_token0').is_in([
                'amg', 'auxiliary_metabolic', 'host_takeover', 'moron', 'anti-defense',
            ])
            | pl.col('empathi_token1').is_in([
                'amg', 'auxiliary_metabolic', 'host_takeover', 'moron', 'anti-defense',
            ])
        ).alias('is_amg_host_takeover'),
    ]).with_columns(
        (
            # Pharokka's "other" is a broad category, not a functional call;
            # require its gene annotation to also be uninformative.
            pl.col('pharokka_category')
                .str.strip_chars()
                .str.to_lowercase()
                .is_in(['', 'unknown function', 'other'])
            & pl.col('pharokka_annot').str.contains(UNINFORMATIVE_ANNOTATION_RE)
            # Phold and Empathi fields can themselves contain functional calls.
            & pl.col('phold_category').str.contains(UNINFORMATIVE_ANNOTATION_RE)
            & pl.col('empathi_annot').str.contains(UNINFORMATIVE_ANNOTATION_RE)
        ).alias('is_unannotated')
    ).with_columns(
        (
            ~pl.col('is_unannotated')
            & ~(
                pl.col('is_capsid_packaging')
                | pl.col('is_dna_metabolism')
                | pl.col('is_tail')
                | pl.col('is_lysis')
                | pl.col('is_connector')
                | pl.col('is_transcription')
                | pl.col('is_integration')
                | pl.col('is_amg_host_takeover')
            )
        ).alias('is_other')
    )


### Defining a breadth of coverage threshold
# Uses the same biological-group membership as parse / discrimination.

gene_breadth_files = glob.glob(
    '../figure_4/activity_profiling/liang_enriched_results/2026-04-03_outputs/'
    'referenceanalyze/coverm/*/bam/*.gene_coverage.tsv.gz'
)

gene_breadth_all = with_bio_group_flags(
    pl.concat([
        pl.read_csv(f, separator='\t').select([
            'breadth', 'pharokka_category', 'pharokka_annot',
            'phold_category', 'empathi_annot',
        ])
        for f in gene_breadth_files
    ])
)

thresholds = np.round(np.arange(0.0, 1.01, 0.05), 2)
CHOSEN_THRESHOLD = GENE_BREADTH_THRESHOLD
PLOT_BIO_GROUP_COLS = BIO_GROUP_COLS + ['unannotated']

rows = []
overall_n = []
for t in thresholds:
    passing = gene_breadth_all.filter(pl.col('breadth') > t)
    n_pass = passing.height
    overall_n.append(n_pass)
    for g in PLOT_BIO_GROUP_COLS:
        prop = None if n_pass == 0 else float(passing[f'is_{g}'].mean())
        rows.append({
            'min_gene_breadth': t,
            'group': g,
            'prop_overall': prop,
            'n_passing': n_pass,
        })

cat_curve_df = pl.DataFrame(rows)
curve_df = pl.DataFrame({
    'min_gene_breadth': thresholds,
    'n_genes_passing': overall_n,
})
print(curve_df)
print(f"\n=== prop_overall by biological group at breadth > {CHOSEN_THRESHOLD} ===")
print(
    cat_curve_df
    .filter(pl.col('min_gene_breadth') == CHOSEN_THRESHOLD)
    .select(['group', 'prop_overall', 'n_passing'])
    .sort('prop_overall', descending=True)
)

fig, ax = plt.subplots(figsize=(8, 5))
for g in PLOT_BIO_GROUP_COLS:
    sub = cat_curve_df.filter(pl.col('group') == g).sort('min_gene_breadth')
    ax.plot(
        sub['min_gene_breadth'].to_list(),
        sub['prop_overall'].to_list(),
        marker='o',
        markersize=3,
        label=g,
        linewidth=2.0 if g != 'unannotated' else 1.0,
        linestyle='--' if g == 'unannotated' else '-',
        color='grey' if g == 'unannotated' else None,
    )
ax.axvline(CHOSEN_THRESHOLD, color='red', linestyle='--', label=f'Chosen threshold ({CHOSEN_THRESHOLD})')
ax.set_xlabel('Minimum gene breadth')
ax.set_ylabel('Proportion of passing genes in group')
ax.set_title('Biological-group composition among genes with breadth > threshold')
ax.set_ylim(0, None)
ax.set_xlim(0, 1)
ax.legend(fontsize=8, loc='best')
plt.tight_layout()
plt.show()


# ===== notebook cell 6 =====
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


# ===== notebook cell 7 =====
# load final metadata
uhvdb_final_metadata = pl.read_csv('../figure_s18/uhvdb_v5_final_metadata_v2.tsv.gz', separator='\t')


# ===== notebook cell 8 =====
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


# ===== notebook cell 9 =====
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


# ===== notebook cell 47 =====
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


# ===== dump training matrix =====
OUT = Path("liang_train_enrich_v_unenrich_final.tsv")
enrich_v_unenrich_final.write_csv(OUT, separator="\t")
print("wrote", OUT, "shape", enrich_v_unenrich_final.shape)
print(enrich_v_unenrich_final["pr_cat"].value_counts())
