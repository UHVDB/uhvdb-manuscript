#!/usr/bin/env python
"""Build UHVDB r4 metadata TSV matching the r5 schema.

Adapted from toolkit/bin/uhvdb_build_metadata.py for a from-scratch build
using packaged figure_1 release tables.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import taxopy

FIG1 = Path("/mmfs1/gscratch/pedslabs_hoffman/carsonjm/CFPhageome/repos/UHVDB/uhvdb-manuscript/figure_1")
FIG3 = Path("/mmfs1/gscratch/pedslabs_hoffman/carsonjm/CFPhageome/repos/UHVDB/uhvdb-manuscript/figure_3")
RELEASE = FIG1 / "uhvdb_human_metag_results" / "uhvdb_2026-03-26-2"
SOURCE_META = FIG1 / "uhvdb_final_metadata.tsv"
LIFESTYLE = FIG3 / "uhvdb_v4_lifestyle.tsv"
CLUSTER_MAPPING = FIG3 / "uhvdb_cluster_mapping.tsv"
R5_METADATA = Path(
    "/mmfs1/gscratch/pedslabs_hoffman/carsonjm/CFPhageome/repos/UHVDB/toolkit2/databases/uhvdb/5.0/uhvdb_metadata.tsv.gz"
)
OUT_PATH = FIG3 / "uhvdb_r4_metadata.tsv.gz"

DEFAULT_TAXDUMP_URL = (
    "https://github.com/shenwei356/gtdb-taxdump/releases/download/v0.6.0/gtdb-taxdump-R226.tar.gz"
)

R5_COLUMNS = [
    "seq_name", "source_db", "body_site", "db_type", "checkv_quality",
    "hash", "seqhash_rep", "uhvdb_id", "genomovar_rep", "genomovar_cluster_id",
    "species_rep", "species_cluster_id", "family_cluster_id", "subfamily_cluster_id",
    "genus_cluster_id", "subgenus_cluster_id", "topology", "n_genes", "genetic_code",
    "virus_score", "n_hallmarks", "contig_length", "provirus", "proviral_length",
    "viral_genes", "host_genes", "completeness", "completeness_method", "kmer_freq",
    "Prediction", "Score", "uhvdb_virus_classification", "aai_expected_length",
    "aai_completeness", "aai_confidence", "aai_error", "aai_num_hits", "aai_top_hit",
    "aai_id", "aai_af", "virus_hallmarks", "plasmid_hallmarks", "genomad_taxonomy",
    "ictv_ref", "ictv_proteinsimilarity", "ictv_class", "ictv_order", "ictv_family",
    "ictv_genus", "ictv_species", "crispr_gtdb_r220_species", "crispr_species_connections",
    "crispr_species_agreement", "crispr_gtdb_r220_genus", "crispr_genus_connections",
    "crispr_genus_agreement", "crispr_gtdb_r220_family", "crispr_family_connections",
    "crispr_family_agreement", "phist_gtdb_r226_species", "phist_species_connections",
    "phist_species_agreement", "phist_gtdb_r226_genus", "phist_genus_connections",
    "phist_genus_agreement", "phist_gtdb_r226_family", "phist_family_connections",
    "phist_family_agreement", "final_host_pred", "host_lineage", "virulent", "temperate",
    "phrog_integrases", "phrog_integration_excision", "empathi_integration",
    "num_proteins", "num_uniprot_ips", "num_card", "num_vfdb", "num_tail", "num_capsid",
    "num_lysis", "mcp_hallmark", "terl_hallmark", "portal_hallmark", "integration_status",
]

ANNOTATION_COMPARE_COLS = [
    "source_db", "body_site", "db_type", "checkv_quality", "hash", "seqhash_rep",
    "topology", "n_genes", "genetic_code", "virus_score", "n_hallmarks", "contig_length",
    "provirus", "proviral_length", "viral_genes", "host_genes", "completeness",
    "completeness_method", "kmer_freq", "Prediction", "Score", "uhvdb_virus_classification",
    "aai_expected_length", "aai_completeness", "aai_confidence", "aai_error", "aai_num_hits",
    "aai_top_hit", "aai_id", "aai_af", "virus_hallmarks", "plasmid_hallmarks",
    "genomad_taxonomy", "ictv_ref", "ictv_proteinsimilarity", "ictv_class", "ictv_order",
    "ictv_family", "ictv_genus", "ictv_species", "crispr_gtdb_r220_species",
    "crispr_species_connections", "crispr_species_agreement", "crispr_gtdb_r220_genus",
    "crispr_genus_connections", "crispr_genus_agreement", "crispr_gtdb_r220_family",
    "crispr_family_connections", "crispr_family_agreement", "phist_gtdb_r226_species",
    "phist_species_connections", "phist_species_agreement", "phist_gtdb_r226_genus",
    "phist_genus_connections", "phist_genus_agreement", "phist_gtdb_r226_family",
    "phist_family_connections", "phist_family_agreement", "final_host_pred", "host_lineage",
    "virulent", "temperate", "phrog_integrases", "phrog_integration_excision",
    "empathi_integration", "num_proteins", "num_uniprot_ips", "num_card", "num_vfdb",
    "num_tail", "num_capsid", "num_lysis", "mcp_hallmark", "terl_hallmark",
    "portal_hallmark", "integration_status",
]

CLUSTER_COMPARE_COLS = [
    "genomovar_rep", "genomovar_cluster_id", "species_rep", "species_cluster_id",
    "family_cluster_id", "subfamily_cluster_id", "genus_cluster_id", "subgenus_cluster_id",
]


def create_base_metadata() -> pl.DataFrame:
    print("Loading source metadata and clustering tables...")
    source_meta = pl.read_csv(SOURCE_META, separator="\t")
    seqhasher = pl.read_csv(RELEASE / "uhvdb_seqhasher.tsv.gz", separator="\t").unique("original_id")
    mapping = pl.read_csv(RELEASE / "uhvdb_mapping.tsv.gz", separator="\t").unique("original_id")
    genomovars = pl.read_csv(RELEASE / "uhvdb_genomovars_info.tsv.gz", separator="\t")
    species = pl.read_csv(RELEASE / "uhvdb_species_info.tsv.gz", separator="\t")
    aaicluster = pl.read_csv(RELEASE / "uhvdb_aaicluster.tsv.gz", separator="\t")

    classify_cols = [
        "uhvdb_id", "topology", "n_genes", "genetic_code", "virus_score",
        "n_hallmarks", "contig_length", "provirus", "proviral_length", "viral_genes",
        "host_genes", "completeness", "completeness_method", "kmer_freq", "Prediction",
        "Score", "uhvdb_virus_classification",
    ]
    classify_df = pl.read_csv(
        RELEASE / "uhvdb_classify.tsv.gz",
        separator="\t",
        columns=classify_cols,
        null_values=["NA", ""],
        infer_schema_length=10000,
    ).unique("uhvdb_id")

    hqfilter_cols = [
        "uhvdb_id", "aai_expected_length", "aai_completeness", "aai_confidence", "aai_error",
        "aai_num_hits", "aai_top_hit", "aai_id", "aai_af",
    ]
    hqfilter_df = pl.read_csv(
        RELEASE / "uhvdb_hqfilter.tsv.gz",
        separator="\t",
        columns=hqfilter_cols,
        null_values=["NA", ""],
        infer_schema_length=10000,
    ).unique("uhvdb_id")

    hcfilter_df = pl.read_csv(
        RELEASE / "uhvdb_hcfilter.tsv.gz",
        separator="\t",
        null_values=["NA"],
        schema_overrides={"virus_hallmarks": pl.Int64, "plasmid_hallmarks": pl.Int64},
    ).unique("uhvdb_id").with_columns([
        pl.col("virus_hallmarks").fill_null(0),
        pl.col("plasmid_hallmarks").fill_null(0),
    ])

    seqhash_reps = (
        seqhasher
        .filter(pl.col("original_id").is_in(set(mapping["original_id"])))
        .unique("original_id")
        .rename({"original_id": "seqhash_rep"})
    )

    # Resolve hash for seq_name, including provirus-trimmed fallback (figure_3a pattern)
    seqhasher_base = seqhasher.with_columns(
        pl.col("original_id")
        .str.extract(r"^(.*)\|provirus", group_index=1)
        .fill_null(pl.col("original_id"))
        .alias("original_id_base")
    )

    print("Building base ID / cluster / classify joins...")
    base = (
        source_meta
        .join(seqhasher.select(["original_id", "hash"]), left_on="seq_name", right_on="original_id", how="left")
        .join(
            seqhasher_base.select(["original_id_base", "hash"]).rename({"hash": "hash_right"}),
            left_on="seq_name",
            right_on="original_id_base",
            how="left",
        )
        .with_columns(
            pl.when(pl.col("hash").is_not_null())
            .then(pl.col("hash"))
            .otherwise(pl.col("hash_right"))
            .alias("hash")
        )
        .drop("hash_right")
        .join(seqhash_reps, on="hash", how="inner")
        .join(mapping, left_on="seqhash_rep", right_on="original_id", how="inner")
        .rename({"new_id": "uhvdb_id"})
        .join(classify_df, on="uhvdb_id", how="left")
        .join(hqfilter_df, on="uhvdb_id", how="left")
        .join(hcfilter_df, on="uhvdb_id", how="left")
        .join(
            genomovars.select(["uhvdb_id", "votu_rep", "cluster_id"]).unique("uhvdb_id"),
            on="uhvdb_id",
            how="left",
        )
        .rename({"votu_rep": "genomovar_rep", "cluster_id": "genomovar_cluster_id"})
        .join(
            species.select(["uhvdb_id", "votu_rep", "cluster_id"]).unique("uhvdb_id"),
            left_on="genomovar_rep",
            right_on="uhvdb_id",
            how="left",
        )
        .rename({"votu_rep": "species_rep", "cluster_id": "species_cluster_id"})
        .join(aaicluster.unique("uhvdb_id"), left_on="species_rep", right_on="uhvdb_id", how="left")
        .unique("seq_name")
        .with_columns([
            pl.col("Score").fill_null(0.0),
            pl.col("completeness").fill_null(0.0),
            pl.col("aai_error").fill_null(0.0),
        ])
    )

    print(
        f"  base rows={base.height:,}  uhvdb_ids={base['uhvdb_id'].n_unique():,}  "
        f"genomovar_reps={(base['uhvdb_id'] == base['genomovar_rep']).sum():,}"
    )

    if CLUSTER_MAPPING.exists():
        cm = pl.read_csv(CLUSTER_MAPPING, separator="\t")
        print(f"  cluster_mapping rows={cm.height:,} (cross-check)")

    return base


def add_taxonomy(metadata_df: pl.DataFrame) -> pl.DataFrame:
    print("Adding ICTV taxonomy...")
    taxonomy_df = (
        pl.read_csv(RELEASE / "uhvdb_ictv_taxonomy.tsv.gz", separator="\t")
        .with_columns(pl.col("normscore").fill_null(0))
        .sort("normscore", descending=True)
        .unique("uhvdb_id", maintain_order=True)
        .rename({
            "taxonomy": "genomad_taxonomy",
            "ref": "ictv_ref",
            "normscore": "ictv_proteinsimilarity",
            "Class": "ictv_class",
            "Order": "ictv_order",
            "Family": "ictv_family",
            "Genus": "ictv_genus",
            "Species": "ictv_species",
        })
    )
    return metadata_df.join(taxonomy_df, left_on="genomovar_rep", right_on="uhvdb_id", how="left")


def add_host_predictions(metadata_df: pl.DataFrame) -> pl.DataFrame:
    print("Adding CRISPR / PHIST host predictions...")
    crisprhost_df = (
        pl.read_csv(RELEASE / "uhvdb_crisprhost.tsv.gz", separator="\t")
        .group_by("uhvdb_id")
        .agg([
            pl.col("top_taxonomy").filter(pl.col("rank") == "species").first().alias("crispr_gtdb_r220_species"),
            pl.col("total_connections").filter(pl.col("rank") == "species").first().alias("crispr_species_connections"),
            pl.col("agreement").filter(pl.col("rank") == "species").first().alias("crispr_species_agreement"),
            pl.col("top_taxonomy").filter(pl.col("rank") == "genus").first().alias("crispr_gtdb_r220_genus"),
            pl.col("total_connections").filter(pl.col("rank") == "genus").first().alias("crispr_genus_connections"),
            pl.col("agreement").filter(pl.col("rank") == "genus").first().alias("crispr_genus_agreement"),
            pl.col("top_taxonomy").filter(pl.col("rank") == "family").first().alias("crispr_gtdb_r220_family"),
            pl.col("total_connections").filter(pl.col("rank") == "family").first().alias("crispr_family_connections"),
            pl.col("agreement").filter(pl.col("rank") == "family").first().alias("crispr_family_agreement"),
        ])
    )

    phisthost_df = (
        pl.read_csv(RELEASE / "uhvdb_phisthost.tsv.gz", separator="\t")
        .group_by("uhvdb_id")
        .agg([
            pl.col("consensus_taxonomy").filter(pl.col("rank") == "species").first().str.replace("s__", "").alias("phist_gtdb_r226_species"),
            pl.col("total_connections").filter(pl.col("rank") == "species").first().alias("phist_species_connections"),
            pl.col("agreement").filter(pl.col("rank") == "species").first().alias("phist_species_agreement"),
            pl.col("consensus_taxonomy").filter(pl.col("rank") == "genus").first().str.replace("g__", "").alias("phist_gtdb_r226_genus"),
            pl.col("total_connections").filter(pl.col("rank") == "genus").first().alias("phist_genus_connections"),
            pl.col("agreement").filter(pl.col("rank") == "genus").first().alias("phist_genus_agreement"),
            pl.col("consensus_taxonomy").filter(pl.col("rank") == "family").first().str.replace("f__", "").alias("phist_gtdb_r226_family"),
            pl.col("total_connections").filter(pl.col("rank") == "family").first().alias("phist_family_connections"),
            pl.col("agreement").filter(pl.col("rank") == "family").first().alias("phist_family_agreement"),
        ])
        .with_columns([
            pl.col("phist_species_agreement").cast(pl.Float64),
            pl.col("phist_genus_agreement").cast(pl.Float64),
            pl.col("phist_family_agreement").cast(pl.Float64),
        ])
    )

    conn_cols = [
        "phist_species_connections", "crispr_species_connections",
        "phist_genus_connections", "crispr_genus_connections",
        "phist_family_connections", "crispr_family_connections",
    ]
    agree_cols = [
        "phist_species_agreement", "crispr_species_agreement",
        "phist_genus_agreement", "crispr_genus_agreement",
        "phist_family_agreement", "crispr_family_agreement",
    ]

    combined_host_df = (
        crisprhost_df
        .join(phisthost_df, on="uhvdb_id", how="full", coalesce=True)
        .with_columns([pl.col(c).cast(pl.Int64) for c in conn_cols])
        .with_columns([pl.col(c).fill_null(0) for c in conn_cols + agree_cols])
        .with_columns(
            pl.when(
                (pl.col("phist_species_connections") >= pl.col("crispr_species_connections"))
                & (pl.col("phist_species_connections") > 0)
            )
            .then(pl.col("phist_gtdb_r226_species"))
            .when(pl.col("phist_species_connections") < pl.col("crispr_species_connections"))
            .then(pl.col("crispr_gtdb_r220_species"))
            .when(
                (pl.col("phist_genus_connections") >= pl.col("crispr_genus_connections"))
                & (pl.col("phist_genus_connections") > 0)
            )
            .then(pl.col("phist_gtdb_r226_genus"))
            .when(pl.col("phist_genus_connections") < pl.col("crispr_genus_connections"))
            .then(pl.col("crispr_gtdb_r220_genus"))
            .when(
                (pl.col("phist_family_connections") >= pl.col("crispr_family_connections"))
                & (pl.col("phist_family_connections") > 0)
            )
            .then(pl.col("phist_gtdb_r226_family"))
            .when(pl.col("phist_family_connections") < pl.col("crispr_family_connections"))
            .then(pl.col("crispr_gtdb_r220_family"))
            .otherwise(None)
            .alias("final_host_pred")
        )
    )

    print("Resolving host lineages with taxopy (GTDB R226)...")
    taxdb_r226 = taxopy.TaxDb(taxdump_url=DEFAULT_TAXDUMP_URL)
    unique_final = (
        combined_host_df.select("final_host_pred").drop_nulls().unique().to_series().to_list()
    )
    rank_prefix = [
        ("superkingdom", "s__"),
        ("phylum", "p__"),
        ("class", "c__"),
        ("order", "o__"),
        ("family", "f__"),
        ("genus", "g__"),
        ("species", "s__"),
    ]
    species_to_lineage: dict[str, str] = {}
    for sp in unique_final:
        taxids = taxopy.taxid_from_name(sp, taxdb_r226)
        if not taxids:
            continue
        taxon = taxopy.Taxon(taxids[0], taxdb_r226)
        rank_name = taxon.rank_name_dictionary
        lineage_parts = [
            f"{prefix}{rank_name[rank]}"
            for rank, prefix in rank_prefix
            if rank in rank_name and rank_name[rank] is not None
        ]
        species_to_lineage[sp] = ";".join(lineage_parts)

    combined_host_df = combined_host_df.with_columns(
        pl.col("final_host_pred")
        .replace(species_to_lineage)
        .str.replace(r"^s__", "d__")
        .alias("host_lineage")
    )

    return metadata_df.join(
        combined_host_df, left_on="genomovar_rep", right_on="uhvdb_id", how="left"
    ).with_columns([
        pl.col(c).fill_null(0)
        for c in conn_cols + agree_cols
    ])


def load_integrated_seqs() -> set[str]:
    print("Loading original integrated-sequence evidence (figure_3d)...")
    uhgv_proviruses = set(
        pl.read_csv(FIG3 / "genomad_viral_stats.tsv", separator="\t", columns=["contig_id", "provirus"])
        .filter(pl.col("provirus") == "Yes")["contig_id"]
    )
    opd_proviruses = set(
        pl.read_csv(FIG3 / "opd_checkv_provirus.tsv", separator="\t", columns=["replace_id"])["replace_id"]
    )
    cgvr_proviruses = set(
        pl.read_csv(FIG3 / "cgvr_checkv_proviruses.tsv", separator="\t", columns=["Viral_bin_id"])["Viral_bin_id"]
    )
    imgvr_proviruses = set(
        pl.read_csv(
            FIG1 / "IMGVR_all_Sequence_information.tsv",
            separator="\t",
            columns=["UVIG", "Topology"],
        )
        .filter(pl.col("Topology") == "Provirus")["UVIG"]
    )
    mmge_proviruses = set(
        pl.read_csv(FIG1 / "all_mge_inf.csv").filter(pl.col("prophage") == True)["MGEs_id"]
    )
    original_classify_proviruses = set(
        pl.read_csv(
            FIG1 / "viruses.csvtk_concat.tsv",
            separator="\t",
            columns=["topology", "provirus", "seq_name"],
        )
        .filter((pl.col("topology") == "Provirus") | (pl.col("provirus") == "Yes"))["seq_name"]
    )
    integrated = (
        uhgv_proviruses
        | opd_proviruses
        | cgvr_proviruses
        | imgvr_proviruses
        | mmge_proviruses
        | original_classify_proviruses
    )
    print(f"  integrated input seqs={len(integrated):,}")
    return integrated


def add_lifestyle_and_integration(metadata_df: pl.DataFrame) -> pl.DataFrame:
    print("Joining lifestyle annotations...")
    lifestyle = pl.read_csv(
        LIFESTYLE,
        separator="\t",
        schema_overrides={
            "phrog_integrases": pl.Int64,
            "phrog_integration_excision": pl.Int64,
            "empathi_integration": pl.Int64,
        },
    ).select([
        "uhvdb_id", "virulent", "temperate",
        "phrog_integrases", "phrog_integration_excision", "empathi_integration",
    ])

    integrated_seqs = load_integrated_seqs()

    # Per-row status (do not propagate member evidence onto genomovar_rep)
    out = (
        metadata_df
        .join(lifestyle, left_on="genomovar_rep", right_on="uhvdb_id", how="left")
        .with_columns(
            pl.when(
                (pl.col("topology") == "Provirus")
                | (pl.col("provirus") == "Yes")
                | pl.col("seq_name").is_in(integrated_seqs)
            )
            .then(pl.lit("integrated"))
            .otherwise(pl.lit("unintegrated"))
            .alias("integration_status")
        )
    )
    return out


def add_protein_aggregates(metadata_df: pl.DataFrame) -> pl.DataFrame:
    print("Computing protein hallmark aggregates (genomovar reps only)...")
    genomovar_reps = set(
        metadata_df.filter(pl.col("uhvdb_id") == pl.col("genomovar_rep"))["uhvdb_id"]
    )
    print(f"  genomovar_reps={len(genomovar_reps):,}")

    prothash_df = (
        pl.read_csv(RELEASE / "uhvdb_proteinhash.tsv.gz", separator="\t")
        .unique("protein_id")
        .with_columns(pl.col("protein_id").str.replace(r"_[^_]*$", "").alias("genomovar_rep"))
        .filter(pl.col("genomovar_rep").is_in(genomovar_reps))
    )
    print(f"  proteinhash rows for reps={prothash_df.height:,}")
    hashes = set(prothash_df["hash"])

    bakta_df = pl.read_csv(
        RELEASE / "uhvdb_bakta.tsv.gz",
        separator="\t",
        columns=["Locus Tag", "Accession"],
        new_columns=["hash", "bakta_acc"],
        null_values=["-"],
        skip_rows=5,
        ignore_errors=True,
    ).filter(pl.col("hash").is_in(hashes))

    foldseek_df = (
        pl.read_csv(
            RELEASE / "uhvdb_foldseek.tsv.gz",
            separator="\t",
            has_header=False,
            columns=["column_1", "column_2"],
            new_columns=["hash", "foldseek_acc"],
        )
        .filter(pl.col("hash").is_in(hashes))
        .unique("hash")
    )

    ips_df = (
        pl.read_csv(
            RELEASE / "uhvdb_interproscan.tsv.gz",
            separator="\t",
            ignore_errors=True,
            has_header=False,
            columns=["column_1", "column_5"],
            new_columns=["hash", "ips_id"],
        )
        .filter(pl.col("hash").is_in(hashes))
        .group_by("hash")
        .agg(pl.col("ips_id").cast(pl.String).str.join(",").alias("ips_id"))
    )

    try:
        card_df = (
            pl.read_csv(
                RELEASE / "uhvdb_card.tsv.gz",
                separator="\t",
                ignore_errors=True,
                has_header=False,
                columns=["column_1", "column_2"],
                new_columns=["hash", "card_acc"],
            )
            .filter(pl.col("hash").is_in(hashes))
            .unique("hash")
        )
    except Exception:
        card_df = pl.DataFrame(schema={"hash": pl.String, "card_acc": pl.String})

    try:
        vfdb_df = (
            pl.read_csv(
                RELEASE / "uhvdb_vfdb.tsv.gz",
                separator="\t",
                ignore_errors=True,
                has_header=False,
                columns=["column_1", "column_2"],
                new_columns=["hash", "vfdb_acc"],
            )
            .filter(pl.col("hash").is_in(hashes))
            .unique("hash")
        )
    except Exception:
        vfdb_df = pl.DataFrame(schema={"hash": pl.String, "vfdb_acc": pl.String})

    pharokka_df = (
        pl.read_csv(
            RELEASE / "uhvdb_pharokka.tsv.gz",
            separator="\t",
            ignore_errors=True,
            columns=["ID", "annot", "category"],
            new_columns=["hash", "pharokka_annot", "pharokka_category"],
        )
        .filter(pl.col("pharokka_annot") != "hypothetical protein")
        .filter(pl.col("hash").is_in(hashes))
        .unique("hash")
    )

    phold_df = (
        pl.read_csv(
            RELEASE / "uhvdb_phold.tsv.gz",
            separator="\t",
            ignore_errors=True,
            columns=["cds_id", "function", "product"],
            new_columns=["hash", "phold_category", "phold_annot"],
        )
        .filter(pl.col("phold_category") != "unknown function")
        .filter(pl.col("hash").is_in(hashes))
        .unique("hash")
    )

    empathi_df = (
        pl.read_csv(
            RELEASE / "uhvdb_empathi.csv.gz",
            ignore_errors=True,
            columns=["", "Annotation"],
            new_columns=["hash", "empathi_annot"],
        )
        .filter(pl.col("hash").is_in(hashes))
        .unique("hash")
    )

    print("  joining protein annotation tables...")
    combined = (
        prothash_df
        .join(bakta_df, on="hash", how="inner")
        .join(foldseek_df, on="hash", how="left")
        .join(ips_df, on="hash", how="left")
        .join(card_df, on="hash", how="left")
        .join(vfdb_df, on="hash", how="left")
        .join(pharokka_df, on="hash", how="left")
        .join(phold_df, on="hash", how="left")
        .join(empathi_df, on="hash", how="left")
    )

    agg = combined.group_by("genomovar_rep").agg([
        pl.len().alias("num_proteins"),
        (
            (pl.col("bakta_acc").is_not_null())
            | (pl.col("foldseek_acc").is_not_null())
            | (pl.col("ips_id").is_not_null())
        ).sum().alias("num_uniprot_ips"),
        pl.col("card_acc").is_not_null().sum().alias("num_card"),
        pl.col("vfdb_acc").is_not_null().sum().alias("num_vfdb"),
        (
            pl.col("pharokka_category").fill_null("").str.contains("tail")
            | pl.col("phold_category").fill_null("").str.contains("tail")
            | pl.col("empathi_annot").fill_null("").str.contains("tail")
        ).sum().alias("num_tail"),
        (
            pl.col("pharokka_category").fill_null("").str.contains("head")
            | pl.col("phold_category").fill_null("").str.contains("head")
            | pl.col("empathi_annot").fill_null("").str.contains("capsid")
        ).sum().alias("num_capsid"),
        (
            pl.col("pharokka_category").fill_null("").str.contains("lysis")
            | pl.col("phold_category").fill_null("").str.contains("lysis")
            | pl.col("empathi_annot").fill_null("").str.contains("lysis")
        ).sum().alias("num_lysis"),
        (
            ((pl.col("pharokka_annot") == "major head protein").sum() == 1)
            | ((pl.col("phold_category") == "major head protein").sum() == 1)
            | ((pl.col("empathi_annot") == "pvp|capsid|major_capsid").sum() == 1)
        ).sum().alias("mcp_hallmark"),
        (
            ((pl.col("pharokka_annot") == "terminase large subunit").sum() == 1)
            | ((pl.col("phold_category") == "terminase large subunit").sum() == 1)
            | ((pl.col("empathi_annot") == "DNA-associated|terminase|packaging_assembly").sum() == 1)
        ).sum().alias("terl_hallmark"),
        (
            ((pl.col("pharokka_annot") == "portal protein").sum() == 1)
            | ((pl.col("phold_category") == "portal protein").sum() == 1)
            | ((pl.col("empathi_annot") == "pvp|portal").sum() == 1)
        ).sum().alias("portal_hallmark"),
    ])

    print(f"  protein aggregates for {agg.height:,} genomovar_reps")
    return metadata_df.join(agg, on="genomovar_rep", how="left")


def finalize_and_write(metadata_df: pl.DataFrame) -> pl.DataFrame:
    import gzip

    print("Finalizing column order and writing output...")
    for col in R5_COLUMNS:
        if col not in metadata_df.columns:
            metadata_df = metadata_df.with_columns(pl.lit(None).alias(col))
    out = metadata_df.select(R5_COLUMNS)

    print(
        f"  rows={out.height:,}  uhvdb_ids={out['uhvdb_id'].n_unique():,}  "
        f"genomovar_reps={out.filter(pl.col('uhvdb_id') == pl.col('genomovar_rep')).height:,}"
    )
    print("  source_db:\n", out["source_db"].value_counts().sort("count", descending=True))
    print("  body_site:\n", out["body_site"].value_counts().sort("count", descending=True))

    with gzip.open(OUT_PATH, "wb") as f:
        out.write_csv(f, separator="\t")
    print(f"Wrote {OUT_PATH}")
    return out


def compare_to_r5(r4: pl.DataFrame) -> None:
    print("Comparing r4 genomovar_rep annotations to r5...")
    # Canonical rows: sequence is the seqhash rep and the UHVDB ID is the genomovar rep
    r4_reps = r4.filter(
        (pl.col("uhvdb_id") == pl.col("genomovar_rep"))
        & (pl.col("seq_name") == pl.col("seqhash_rep"))
    ).unique("uhvdb_id")
    r5 = pl.read_csv(
        R5_METADATA,
        separator="\t",
        null_values=["NA", ""],
        infer_schema_length=10000,
    ).filter(
        (pl.col("uhvdb_id") == pl.col("genomovar_rep"))
        & (pl.col("seq_name") == pl.col("seqhash_rep"))
    ).unique("uhvdb_id")

    joined = r4_reps.join(r5, on="uhvdb_id", how="inner", suffix="_r5")
    print(f"  r4 canonical genomovar_reps={r4_reps.height:,}  overlapping in r5={joined.height:,}")

    fill0_int = [
        "virus_hallmarks", "plasmid_hallmarks",
        "crispr_species_connections", "crispr_genus_connections", "crispr_family_connections",
        "phist_species_connections", "phist_genus_connections", "phist_family_connections",
        "phrog_integrases", "phrog_integration_excision", "empathi_integration",
    ]
    fill0_float = [
        "Score", "completeness", "aai_error",
        "crispr_species_agreement", "crispr_genus_agreement", "crispr_family_agreement",
        "phist_species_agreement", "phist_genus_agreement", "phist_family_agreement",
    ]
    for c in fill0_int:
        if c in joined.columns and f"{c}_r5" in joined.columns:
            joined = joined.with_columns(pl.col(c).fill_null(0), pl.col(f"{c}_r5").fill_null(0))
    for c in fill0_float:
        if c in joined.columns and f"{c}_r5" in joined.columns:
            joined = joined.with_columns(pl.col(c).fill_null(0.0), pl.col(f"{c}_r5").fill_null(0.0))

    print("\n  Annotation column match rates:")
    match_exprs = []
    for col in ANNOTATION_COMPARE_COLS:
        if col not in r4_reps.columns:
            continue
        r5_col = f"{col}_r5"
        if r5_col not in joined.columns:
            continue
        match_expr = pl.col(col).cast(pl.String).eq_missing(pl.col(r5_col).cast(pl.String))
        rate = joined.select(match_expr.mean()).item()
        n_mismatch = joined.select((~match_expr).sum()).item()
        print(f"    {col}: {rate:.4%} match ({n_mismatch:,} mismatches)")
        match_exprs.append(match_expr)

    if match_exprs:
        fully_matching = match_exprs[0]
        for expr in match_exprs[1:]:
            fully_matching = fully_matching & expr
        n_full = joined.select(fully_matching.sum()).item()
        print(
            f"\n  Fully matching annotation rows: {n_full:,} / {joined.height:,} "
            f"({n_full / joined.height:.4%})"
        )

    print("\n  Clustering column match rates (informational; may differ after r5 increment):")
    for col in CLUSTER_COMPARE_COLS:
        r5_col = f"{col}_r5"
        if r5_col not in joined.columns:
            continue
        match_expr = pl.col(col).cast(pl.String).eq_missing(pl.col(r5_col).cast(pl.String))
        rate = joined.select(match_expr.mean()).item()
        print(f"    {col}: {rate:.4%} match")

    for col in ["host_lineage", "integration_status", "num_proteins", "ictv_family", "final_host_pred"]:
        r5_col = f"{col}_r5"
        if r5_col not in joined.columns:
            continue
        mism = joined.filter(
            ~pl.col(col).cast(pl.String).eq_missing(pl.col(r5_col).cast(pl.String))
        ).select(["uhvdb_id", col, r5_col]).head(5)
        if mism.height:
            print(f"\n  Sample mismatches for {col}:")
            print(mism)


def main() -> None:
    base = create_base_metadata()
    with_tax = add_taxonomy(base)
    with_host = add_host_predictions(with_tax)
    with_life = add_lifestyle_and_integration(with_host)
    with_prot = add_protein_aggregates(with_life)
    r4 = finalize_and_write(with_prot)
    compare_to_r5(r4)


if __name__ == "__main__":
    main()
