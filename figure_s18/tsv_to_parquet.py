#!/usr/bin/env python3

import argparse
import polars as pl

def main():
    parser = argparse.ArgumentParser(description="Convert TSV to Parquet with Polars using streaming")
    parser.add_argument("input_tsv", help="Input TSV file")
    parser.add_argument("output_parquet", help="Output Parquet file")
    args = parser.parse_args()

    lf = pl.scan_csv(
        args.input_tsv,
        separator="\t",
        infer_schema_length=1000,
        ignore_errors=True,
    )

    lf.sink_parquet(args.output_parquet)

if __name__ == "__main__":
    main()