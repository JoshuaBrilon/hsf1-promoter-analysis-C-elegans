from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


def main():
    root = Path(__file__).resolve().parents[1]

    genome_fasta = root / "data" / "raw" / "libuda_N2_genome.fasta"
    input_path = root / "data" / "processed" / "promoters_with_gene_names.txt"

    dedup_txt = root / "data" / "processed" / "promoters_with_gene_names_dedup.txt"
    positive_fasta = root / "data" / "processed" / "positivePromoters_dedup.fasta"
    negative_fasta = root / "data" / "processed" / "negativePromoters_dedup.fasta"
    combined_fasta = root / "data" / "processed" / "allPromoters_dedup.fasta"

    print("Loading promoters...", flush=True)
    df = pd.read_csv(input_path, sep="\t", header=None)

    original_count = len(df)
    df_dedup = df.drop_duplicates()
    dedup_count = len(df_dedup)

    df_dedup.to_csv(dedup_txt, sep="\t", index=False, header=False)

    print(f"Original rows: {original_count}", flush=True)
    print(f"After deduplication: {dedup_count}", flush=True)
    print(f"Removed duplicates: {original_count - dedup_count}", flush=True)
    print(f"Saved: {dedup_txt}", flush=True)

    print("Loading genome...", flush=True)
    genome = {}
    for record in SeqIO.parse(genome_fasta, "fasta"):
        genome[record.id] = str(record.seq)

    print("Building positive-strand promoter FASTA...", flush=True)
    with open(dedup_txt, "r") as infile, open(positive_fasta, "w") as outfile:
        for line in infile:
            columns = line.strip().split("\t")

            chrom = columns[0]
            start = int(columns[1])
            end = int(columns[2])
            gene_name = columns[3]
            strand = columns[5]

            if strand == "+":
                seq = genome[chrom][start:end]
                outfile.write(f">{gene_name}|{chrom}:{start}-{end}({strand})\n")
                outfile.write(f"{seq}\n")

    print("Building negative-strand promoter FASTA...", flush=True)
    with open(dedup_txt, "r") as infile, open(negative_fasta, "w") as outfile:
        for line in infile:
            columns = line.strip().split("\t")

            chrom = columns[0]
            start = int(columns[1])
            end = int(columns[2])
            gene_name = columns[3]
            strand = columns[5]

            if strand == "-":
                seq = genome[chrom][start:end]
                seq = str(Seq(seq).reverse_complement())
                outfile.write(f">{gene_name}|{chrom}:{start}-{end}({strand})\n")
                outfile.write(f"{seq}\n")

    print("Combining FASTA files...", flush=True)
    with open(combined_fasta, "w") as outfile:
        with open(positive_fasta, "r") as infile:
            outfile.write(infile.read())
        with open(negative_fasta, "r") as infile:
            outfile.write(infile.read())

    print(f"Saved: {positive_fasta}", flush=True)
    print(f"Saved: {negative_fasta}", flush=True)
    print(f"Saved: {combined_fasta}", flush=True)
    print()
    print("Next commands:", flush=True)
    print(
        "fasta-get-markov -m 1 data/processed/allPromoters_dedup.fasta data/processed/background_dedup.txt",
        flush=True,
    )
    print(
        "fimo --oc fimo_out_dedupes --bgfile data/processed/background_dedup.txt data/motif/MA2169.1.meme data/processed/allPromoters_dedup.fasta",
        flush=True,
    )

if __name__ == "__main__":
    main()