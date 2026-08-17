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
