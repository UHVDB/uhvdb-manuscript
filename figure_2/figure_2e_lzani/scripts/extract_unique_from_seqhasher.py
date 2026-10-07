#!/usr/bin/env python3
"""Pick one sequence ID per unique hash from seqhasher TSV chunks.

Writes a plain ID list for seqkit grep (does not touch sequence columns).
Seqhasher TSV columns (no header): seq_id, hash, sequence
"""

import argparse
import gzip
import sys
from pathlib import Path


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seqhasher-dir", required=True, type=Path)
    parser.add_argument("--output-ids", required=True, type=Path)
    parser.add_argument("--glob", default="*.seqhasher.tsv.gz")
    args = parser.parse_args()

    chunks = sorted(args.seqhasher_dir.glob(args.glob))
    if not chunks:
        sys.exit("No seqhasher chunks matching %s in %s" % (args.glob, args.seqhasher_dir))

    args.output_ids.parent.mkdir(parents=True, exist_ok=True)
    seen = set()
    n_in = 0
    n_out = 0

    with open(args.output_ids, "w") as out:
        for chunk in chunks:
            print("reading %s" % chunk.name, flush=True)
            with open_text(chunk) as fh:
                for line in fh:
                    # Only need id + hash; avoid splitting the huge sequence field
                    tab1 = line.find("\t")
                    if tab1 < 0:
                        continue
                    tab2 = line.find("\t", tab1 + 1)
                    if tab2 < 0:
                        continue
                    seq_id = line[:tab1]
                    hsh = line[tab1 + 1 : tab2]
                    n_in += 1
                    if hsh in seen:
                        continue
                    seen.add(hsh)
                    out.write(seq_id + "\n")
                    n_out += 1
                    if n_out % 100000 == 0:
                        print("  unique ids: {:,}".format(n_out), flush=True)

    print("done: rows={:,} unique_ids={:,} -> {}".format(n_in, n_out, args.output_ids), flush=True)


if __name__ == "__main__":
    main()
