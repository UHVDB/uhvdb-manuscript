#!/usr/bin/env python

import argparse
import gzip
import os
import sys


def parse_args(args=None):
    description = "Identify terminal repeats in FASTA sequences."
    epilog = "Example usage: python trfinder.py --input sequences.fasta.gz --prefix test"

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-i",
        "--input",
        help="Path to input FASTA (gzipped) to be searched for terminal repeats.",
    )
    parser.add_argument(
        "--min_dtr", help="Minimum DTR length to be considered a direct terminal repeat.", default=20, type=int
    )
    parser.add_argument(
        "--min_itr", help="Minimum ITR length to be considered an inverse terminal repeat.", default=20, type=int
    )
    parser.add_argument(
        "--max_itr", help="Maximum ITR length to be considered an inverse terminal repeat.", default=1000, type=int
    )
    parser.add_argument(
        "--trim_dtrs", help="Trim DTRs from the end of DTR contigs.", action='store_true'
    )
    parser.add_argument(
        "-p",
        "--prefix",
        help="Prefix to be used in naming output files.",
    )
    parser.add_argument(
        "-s",
        "--source_db",
        help="Abbreviation for the database the contigs were derived from.",
    )
    parser.add_argument(
        "-w",
        "--write_all",
        help="Write all contigs regardless of TR status.",
        type=str
    )
    parser.add_argument(
        "-n",
        "--use_contig_name",
        help="Use original contig name when generating TR IDs.",
    )
    return parser.parse_args(args)


def parse_fasta(handle):
    try:
        id = next(handle).strip()[1:]
    except:
        return
    seq = ""
    for line in handle:
        if line[0] == ">":
            yield id, seq
            id = line.strip()[1:]
            seq = ""
        else:
            seq += line.rstrip().upper()
    yield id, seq


def fetch_dtr(fullseq, min_length):
    startseq = fullseq[:min_length]
    # Find the first occurrence of startseq in the second half of fullseq
    startpos = fullseq.find(startseq, len(fullseq) // 2)
    while startpos != -1:
        endseq = fullseq[startpos:]
        # Check if the match extends to the contig end
        if fullseq.startswith(endseq):
            return endseq
        # Find the next occurrence of startseq in the second half
        startpos = fullseq.find(startseq, startpos + 1)
    return None


def reverse_complement(seq):
    trans = str.maketrans("ACTG", "TGAC")
    return seq[::-1].translate(trans)


def fetch_itr(seq, min_len, max_len):
    rev = reverse_complement(seq)
    # see if minimal substring occurs at end
    if seq[:min_len] == rev[:min_len]:
        # extend to maximum substring, up to <max_len>
        i = min_len + 1
        while seq[:i] == rev[:i] and i <= max_len:
            i += 1
        return seq[: i - 1]
    # no match
    else:
        return None


def fetch_trs(input, source_db, prefix, min_dtr, min_itr, max_itr, write_all, use_contig_name):
    trs = []
    try:
        if not os.path.exists(input):
            return trs
        infile = gzip.open(input, "rt") if input.split(".")[-1] == "gz" else open(input)
        handle = parse_fasta(infile)
        for contig_index, contig in enumerate(handle):
            contig_id, contig_seq = contig
            dtr = fetch_dtr(contig_seq, min_dtr)
            itr = fetch_itr(contig_seq, min_itr, max_itr)
            if not write_all and not dtr and not itr:
                continue
            tr = {
                "sample_id": source_db + prefix,
                "contig_id": source_db + prefix + "_" + str(contig_index + 1),
                "contig_name": contig_id,
                "contig_len": len(contig_seq),
            }
            if dtr:
                tr.update({"tr_type": "DTR", "tr_seq": dtr})
                tr_nt_counts = [tr["tr_seq"].count(_) for _ in list("ACGT")]
                tr["contig_seq"] = contig_seq
                tr["trimmed_contig_seq"] = contig_seq[:-len(dtr)]
                tr["trimmed_contig_len"] = len(contig_seq[:-len(dtr)])
                tr["tr_len"] = len(tr["tr_seq"])
                tr["tr_nt_acgt_count"] = sum(tr_nt_counts)
                tr["tr_nt_max_freq"] = (
                    round(100.0 * max(tr_nt_counts) / sum(tr_nt_counts), 2) if sum(tr_nt_counts) > 0 else None
                )
                tr["tr_nt_n_count"] = tr["tr_len"] - sum(tr_nt_counts)
                trs.append(tr)
            elif itr:
                tr.update({"tr_type": "ITR", "tr_seq": itr})
                tr_nt_counts = [tr["tr_seq"].count(_) for _ in list("ACGT")]
                tr["contig_seq"] = contig_seq
                tr["trimmed_contig_seq"] = contig_seq
                tr["trimmed_contig_len"] = len(contig_seq)
                tr["tr_len"] = len(tr["tr_seq"])
                tr["tr_nt_acgt_count"] = sum(tr_nt_counts)
                tr["tr_nt_max_freq"] = (
                    round(100.0 * max(tr_nt_counts) / sum(tr_nt_counts), 2) if sum(tr_nt_counts) > 0 else None
                )
                tr["tr_nt_n_count"] = tr["tr_len"] - sum(tr_nt_counts)
                trs.append(tr)
            else:
                tr.update({"tr_type": "None", "tr_seq": None})
                tr["contig_seq"] = contig_seq
                tr_nt_counts = 0
                tr["tr_len"] = 0
                tr["tr_nt_acgt_count"] = 0
                tr["tr_nt_max_freq"] = 0
                tr["tr_nt_n_count"] = 0
                trs.append(tr)
            if use_contig_name:
                tr["contig_id"] = tr["contig_name"]
    except:
        pass
    return trs


def write_batch(sample, trim_dtrs, results):
    with open(f"./{sample}.trfinder.tsv", "w") as tsv_outfile:
        with open(f"./{sample}.trfinder.fna", "w") as fasta_outfile:
            fields = [
                "contig_id",
                "contig_name",
                # 'contig_seq',
                "contig_len",
                "tr_type",
                "tr_seq",
                "tr_len",
                "tr_nt_acgt_count",
                "tr_nt_n_count",
                "tr_nt_max_freq",
            ]
            if trim_dtrs:
                fields.append("trimmed_contig_len")
            tsv_outfile.write("\t".join(fields) + "\n")
            for _ in results:
                # write TR info to tsv
                tsv_outfile.write("\t".join([str(_[f]) for f in fields]) + "\n")
                # write TR sequence to fasta
                if trim_dtrs:
                    fasta_outfile.write(">" + str(_["contig_id"]) + "\n" + str(_["trimmed_contig_seq"]) + "\n")
                else:
                    fasta_outfile.write(">" + str(_["contig_id"]) + "\n" + str(_["contig_seq"]) + "\n")


def main(args=None):
    args = parse_args(args)
    results = fetch_trs(
        args.input,
        args.source_db,
        args.prefix,
        args.min_dtr,
        args.min_itr,
        args.max_itr,
        args.write_all,
        args.use_contig_name)

    if args.trim_dtrs:
        trim_dtrs = True
    else:
        trim_dtrs = False

    write_batch(args.prefix, trim_dtrs, results)


if __name__ == "__main__":
    sys.exit(main())
