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
