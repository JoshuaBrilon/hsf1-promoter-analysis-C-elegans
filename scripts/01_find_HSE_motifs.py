from pathlib import Path
from Bio import SeqIO
import re
import matplotlib.pyplot as plt
import os

PATTERN_FORWARD = r"GAA..TTC..GAA"
PATTERN_REVERSE = r"TTC..GAA..TTC"

COLORS = ["blue", "green", "red", "orange", "purple", "cyan"]
ROMAN_CHROMS = ["I", "II", "III", "IV", "V", "X"]


def find_hse_motifs(sequence: str):
    positions = []
    hits = []

    forward_count = 0
    reverse_count = 0

    for match in re.finditer(PATTERN_FORWARD, sequence):
        start = match.start()
        end = match.end() - 1
        positions.append(start)
        hits.append((start, end, "+", match.group()))
        forward_count += 1

    for match in re.finditer(PATTERN_REVERSE, sequence):
        start = match.start()
        end = match.end() - 1
        positions.append(start)
        hits.append((start, end, "-", match.group()))
        reverse_count += 1

    hits.sort(key=lambda x: x[0])
    positions.sort()

    return positions, hits, forward_count, reverse_count


def save_hits(file_path: Path, chrom_id: str, hits):
    with open(file_path, "w") as f:
        for start, end, strand, motif in hits:
            f.write(f"{chrom_id}\t{start}\t{end}\t{strand}\t{motif}\n")


def plot_histogram(positions, color, title, output_path):
    plt.hist(positions, bins=200, color=color)
    plt.title(title)
    plt.xlabel("Genomic Position")
    plt.ylabel("Frequency")
    plt.savefig(output_path)
    plt.clf()


def main():
    ROOT = Path(__file__).resolve().parents[1]

    fasta_path = ROOT / "data" / "raw" / "libuda_N2_genome.fasta"

    output_folder = ROOT / "data" / "processed" / "chromosomes_01"
    plot_folder = ROOT / "results" / "plots"
    bounds_file = ROOT / "data" / "processed" / "chrom_bounds.txt"
    all_hses_path = ROOT / "data" / "processed" / "allHSEs.txt"

    output_folder.mkdir(parents=True, exist_ok=True)
    plot_folder.mkdir(parents=True, exist_ok=True)

    all_positions_by_chrom = []
    chrom_bounds = []

    with open(all_hses_path, "w") as all_hses_file:
        for i, chromosome in enumerate(SeqIO.parse(fasta_path, "fasta")):
            chrom_id = chromosome.id.strip()
            sequence = str(chromosome.seq)

            positions, hits, fwd, rev = find_hse_motifs(sequence)

            print(f"{chrom_id}: forward={fwd}, reverse={rev}, total={fwd + rev}")

            save_hits(output_folder / f"{chrom_id}.txt", chrom_id, hits)

            for start, end, strand, motif in hits:
                all_hses_file.write(f"{chrom_id}\t{start}\t{end}\t{strand}\t{motif}\n")

            plot_histogram(
                positions,
                COLORS[i],
                f"Chromosome {chrom_id}",
                plot_folder / f"{chrom_id}_HSF_histogram.png"
            )

            all_positions_by_chrom.append(positions)
            chrom_bounds.append((chrom_id, len(chromosome.seq)))

    for i, positions in enumerate(all_positions_by_chrom):
        plt.hist(
            positions,
            bins=200,
            color=COLORS[i],
            alpha=0.5,
            label=f"Chr {ROMAN_CHROMS[i]}"
        )

    plt.legend()
    plt.title("HSF Binding Across Chromosomes")
    plt.xlabel("Genomic Position")
    plt.ylabel("Frequency")
    plt.savefig(plot_folder / "HSF_Binding_Total.png")
    plt.clf()

    with open(bounds_file, "w") as f:
        for chrom_id, length in chrom_bounds:
            f.write(f"{chrom_id}\t1\t{length}\n")


if __name__ == "__main__":
    main()