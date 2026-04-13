from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq


def main() -> None:
    root = Path(__file__).resolve().parents[1]

    genome_fasta = root / "data" / "raw" / "libuda_N2_genome.fasta"
    promoter_file = root / "data" / "processed" / "promoters_with_gene_names.txt"

    positive_output = root / "data" / "processed" / "positivePromoters.fasta"
    negative_output = root / "data" / "processed" / "negativePromoters.fasta"
    combined_output = root / "data" / "processed" / "allPromoters.fasta"

    genome = {}
    for record in SeqIO.parse(genome_fasta, "fasta"):
        genome[record.id] = str(record.seq)

    # positive strand promoters
    with open(promoter_file, "r") as infile, open(positive_output, "w") as outfile:
        for line in infile:
            columns = line.strip().split("\t")

            if columns[5] == "+":
                chrom = columns[0]
                start = int(columns[1])
                end = int(columns[2])
                strand = columns[5]

                cut_seq = genome[chrom][start:end]

                outfile.write(f">{chrom}:{start}-{end}({strand})\n")
                outfile.write(f"{cut_seq}\n")

    # negative strand promoters
    with open(promoter_file, "r") as infile, open(negative_output, "w") as outfile:
        for line in infile:
            columns = line.strip().split("\t")

            if columns[5] == "-":
                chrom = columns[0]
                start = int(columns[1])
                end = int(columns[2])
                strand = columns[5]

                seq = genome[chrom][start:end]
                seq = str(Seq(seq).reverse_complement())

                outfile.write(f">{chrom}:{start}-{end}({strand})\n")
                outfile.write(f"{seq}\n")

    # combine both files
    with open(combined_output, "w") as outfile:
        with open(positive_output, "r") as infile:
            outfile.write(infile.read())

        with open(negative_output, "r") as infile:
            outfile.write(infile.read())

    print(f"Saved: {positive_output}", flush=True)
    print(f"Saved: {negative_output}", flush=True)
    print(f"Saved: {combined_output}", flush=True)
    print()
    print("Example FIMO command:", flush=True)
    print(
        "fimo --bgfile background.txt data/motif/MA2169.1.meme data/processed/allPromoters.fasta",
        flush=True,
    )


if __name__ == "__main__":
    main()