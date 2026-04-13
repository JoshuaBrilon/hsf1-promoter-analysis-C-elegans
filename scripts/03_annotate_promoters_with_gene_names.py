from pathlib import Path

import pandas as pd


PROMOTER_LENGTH = 1000


def load_gff(gff_path: Path) -> pd.DataFrame:
    return pd.read_csv(
        gff_path,
        sep="\t",
        header=None,
        names=[
            "chromosome",
            "program",
            "feature",
            "start",
            "stop",
            "score",
            "strand",
            "phase",
            "info",
        ],
        comment="#",
    )


def load_gene_ids(gene_ids_path: Path) -> pd.DataFrame:
    gene_id = pd.read_csv(
        gene_ids_path,
        sep=",",
        comment="#",
        header=None,
        names=[
            "taxon_id",
            "WBGene_ID",
            "gene_name",
            "clone_name",
            "living_status",
            "coding_status",
        ],
    )
    return gene_id[["WBGene_ID", "gene_name"]]


def build_promoters_with_gene_ids(gff: pd.DataFrame) -> pd.DataFrame:
    mrna = gff[gff["feature"] == "mRNA"].copy()

    plus_strand = mrna[mrna["strand"] == "+"].copy()
    minus_strand = mrna[mrna["strand"] == "-"].copy()

    plus_strand["promoter_start"] = (plus_strand["start"] - PROMOTER_LENGTH).clip(lower=0)
    plus_strand["promoter_stop"] = plus_strand["start"]

    minus_strand["promoter_start"] = minus_strand["stop"]
    minus_strand["promoter_stop"] = minus_strand["stop"] + PROMOTER_LENGTH

    full_strands = pd.concat([plus_strand, minus_strand], ignore_index=True)

    promoters = full_strands[
        ["chromosome", "promoter_start", "promoter_stop", "strand", "info"]
    ].copy()

    promoters["WBGene_ID"] = promoters["info"].str.extract(r"Parent=gene:(WBGene\d+)")

    return promoters


def annotate_promoters(promoters: pd.DataFrame, gene_lookup: pd.DataFrame) -> pd.DataFrame:
    annotated = promoters.merge(gene_lookup, on="WBGene_ID", how="left")

    annotated = annotated[
        ["chromosome", "promoter_start", "promoter_stop", "gene_name", "strand"]
    ].copy()

    annotated["score"] = "."

    annotated = annotated[
        ["chromosome", "promoter_start", "promoter_stop", "gene_name", "score", "strand"]
    ].copy()

    return annotated


def main() -> None:
    root = Path(__file__).resolve().parents[1]

    gff_path = root / "data" / "raw" / "N2.genome.annotations.gff3"
    gene_ids_path = root / "data" / "raw" / "C_elegans.current.geneIDs.txt"

    processed_dir = root / "data" / "processed"
    processed_dir.mkdir(parents=True, exist_ok=True)

    output_path = processed_dir / "promoters_with_gene_names.txt"
    all_hses_input = processed_dir / "allHSEs.txt"
    all_hses_output = processed_dir / "allHSEs_bed6.txt"

    print("Loading GFF3...", flush=True)
    gff = load_gff(gff_path)

    print("Loading gene ID table...", flush=True)
    gene_lookup = load_gene_ids(gene_ids_path)

    print("Building promoter regions...", flush=True)
    promoters = build_promoters_with_gene_ids(gff)

    print("Annotating promoters with gene names...", flush=True)
    annotated_promoters = annotate_promoters(promoters, gene_lookup)

    promoters_with_gene_names = annotated_promoters.dropna(subset=["gene_name"]).copy()
    promoters_with_gene_names.to_csv(output_path, sep="\t", index=False, header=False)

    print("Converting allHSEs.txt...", flush=True)
    hses = pd.read_csv(
        all_hses_input,
        sep="\t",
        header=None,
        names=["chromosome", "start", "end", "strand", "motif"],
    )

    # Reorder to BED6-style:
    # chrom, start, end, name, score, strand
    hses["score"] = "."
    hses = hses[["chromosome", "start", "end", "motif", "score", "strand"]]
    hses.to_csv(all_hses_output, sep="\t", index=False, header=False)

    print(f"Total promoter regions: {len(annotated_promoters)}", flush=True)
    print(f"Promoters with matched gene names: {len(promoters_with_gene_names)}", flush=True)
    print(f"Saved: {output_path}", flush=True)
    print(f"Saved: {all_hses_output}", flush=True)
    print()
    print("Run this from the repo root:", flush=True)
    print(
        "bedtools intersect -a data/processed/allHSEs_bed6.txt "
        "-b data/processed/promoters_with_gene_names.txt "
        "-wa -wb -s > data/processed/promoter_overlaps_with_all_genes.txt",
        flush=True,
    )


if __name__ == "__main__":
    main()