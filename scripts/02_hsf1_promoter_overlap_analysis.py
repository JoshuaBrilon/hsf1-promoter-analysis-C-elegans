from pathlib import Path
import re

import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt


PATTERN_FORWARD = re.compile(r"GAA..TTC..GAA")
PATTERN_REVERSE = re.compile(r"TTC..GAA..TTC")
PROMOTER_LENGTH = 1000
COLORS = ["blue", "green", "red", "orange", "purple", "cyan"]


def find_hse_hits(fasta_path: Path, pattern: re.Pattern, strand_label: str) -> pd.DataFrame:
    """Scan a genome FASTA for motif matches and return BED-like coordinates."""
    hits = []

    for record in SeqIO.parse(fasta_path, "fasta"):
        chrom_id = record.id.strip()
        sequence = str(record.seq)

        for match in pattern.finditer(sequence):
            start = match.start()
            end = match.end() - 1
            hits.append([chrom_id, start, end, strand_label])

    return pd.DataFrame(hits, columns=["chromosome", "start", "end", "strand"])


def load_gff(gff_path: Path) -> pd.DataFrame:
    """Load GFF3 annotation file."""
    return pd.read_csv(
        gff_path,
        sep="\t",
        header=None,
        names=["chromosome", "program", "feature", "start", "stop", "score", "strand", "phase", "info"],
        comment="#",
    )


def build_promoters(gff: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build TSS table and strand-specific promoter regions from mRNA annotations."""
    mrna = gff[gff["feature"] == "mRNA"].copy()

    mrna_plus = mrna[mrna["strand"] == "+"].copy()
    mrna_minus = mrna[mrna["strand"] == "-"].copy()

    # Plus strand
    mrna_plus["TSS"] = mrna_plus["start"]
    mrna_plus["TSE"] = mrna_plus["stop"]
    mrna_plus["promoter_start"] = (mrna_plus["start"] - PROMOTER_LENGTH).clip(lower=0)
    mrna_plus["promoter_end"] = mrna_plus["start"]

    # Minus strand
    mrna_minus["TSS"] = mrna_minus["stop"]
    mrna_minus["TSE"] = mrna_minus["start"]
    mrna_minus["promoter_start"] = mrna_minus["stop"]
    mrna_minus["promoter_end"] = mrna_minus["stop"] + PROMOTER_LENGTH

    tss = pd.concat(
        [
            mrna_plus[["chromosome", "TSS", "TSE", "strand"]],
            mrna_minus[["chromosome", "TSS", "TSE", "strand"]],
        ],
        ignore_index=True,
    )

    promoter_plus = mrna_plus[["chromosome", "promoter_start", "promoter_end", "strand"]].copy()
    promoter_minus = mrna_minus[["chromosome", "promoter_start", "promoter_end", "strand"]].copy()

    return tss, promoter_plus, promoter_minus


def plot_overlap_histograms(overlaps_path: Path, output_dir: Path) -> None:
    """Plot chromosome-specific histograms of HSF hits overlapping promoter regions."""
    overlaps = pd.read_csv(overlaps_path, sep="\t", header=None)

    print(f"Overlap file shape: {overlaps.shape}", flush=True)

    # Expected bedtools output:
    # hsf_chr hsf_start hsf_end prom_chr prom_start prom_end prom_strand
    if overlaps.shape[1] == 7:
        overlaps.columns = [
            "hsf_chr",
            "hsf_start",
            "hsf_end",
            "prom_chr",
            "prom_start",
            "prom_end",
            "prom_strand",
        ]
    else:
        raise ValueError(f"Unexpected number of columns in overlap file: {overlaps.shape[1]}")

    overlaps["hsf_start"] = pd.to_numeric(overlaps["hsf_start"], errors="coerce")
    overlaps = overlaps.dropna(subset=["hsf_chr", "hsf_start"])

    output_dir.mkdir(parents=True, exist_ok=True)

    chromosomes = ["N2_chrI", "N2_chrII", "N2_chrIII", "N2_chrIV", "N2_chrV", "N2_chrX"]

    print("Unique values in hsf_chr:")
    print(overlaps["hsf_chr"].unique(), flush=True)

    for i, chrom in enumerate(chromosomes):
        chrom_data = overlaps[overlaps["hsf_chr"] == chrom]
        positions = chrom_data["hsf_start"]

        print(f"Plotting {chrom}: {len(positions)} hits", flush=True)

        if positions.empty:
            continue

        plt.figure()
        plt.hist(positions, bins=100, color=COLORS[i % len(COLORS)])
        plt.title(f"HSF Binding Sites in {chrom} within promoter regions")
        plt.xlabel("Genomic Position")
        plt.ylabel("Count")
        plt.savefig(output_dir / f"{chrom}_HSF_promoter_regions_histogram.png")
        plt.close()

def main() -> None:
    print("Starting script...", flush=True)

    root = Path(__file__).resolve().parents[1]
    print(f"Project root: {root}", flush=True)

    fasta_path = root / "data" / "raw" / "libuda_N2_genome.fasta"
    gff_path = root / "data" / "raw" / "N2.genome.annotations.gff3"
    print(f"FASTA path: {fasta_path}", flush=True)
    print(f"GFF path: {gff_path}", flush=True)

    processed_dir = root / "data" / "processed"
    results_dir = root / "results"

    forward_hits_path = processed_dir / "hsf1_forward_hits.txt"
    reverse_hits_path = processed_dir / "hsf1_reverse_hits.txt"
    tss_path = processed_dir / "mRNA_TSS.txt"
    promoter_plus_path = processed_dir / "mRNApromoterSeqPlus.txt"
    promoter_minus_path = processed_dir / "mRNApromoterSeqMinus.txt"

    overlaps_path = processed_dir / "promoter_overlaps_all.txt"
    histogram_dir = results_dir / "plots" / "promoter_overlap_histograms"

    processed_dir.mkdir(parents=True, exist_ok=True)
    histogram_dir.mkdir(parents=True, exist_ok=True)

    print("Finding forward hits...", flush=True)
    forward_hits = find_hse_hits(fasta_path, PATTERN_FORWARD, "+")
    print("Forward hits finished.", flush=True)

    print("Finding reverse hits...", flush=True)
    reverse_hits = find_hse_hits(fasta_path, PATTERN_REVERSE, "-")
    print("Reverse hits finished.", flush=True)

    forward_hits[["chromosome", "start", "end"]].to_csv(
        forward_hits_path, sep="\t", index=False, header=False
    )
    reverse_hits[["chromosome", "start", "end"]].to_csv(
        reverse_hits_path, sep="\t", index=False, header=False
    )

    print(f"Forward hits: {len(forward_hits)}", flush=True)
    print(f"Reverse hits: {len(reverse_hits)}", flush=True)

    print("Loading GFF...", flush=True)
    gff = load_gff(gff_path)
    print("GFF loaded.", flush=True)

    print("Building promoter regions...", flush=True)
    tss, promoter_plus, promoter_minus = build_promoters(gff)
    print("Promoters built.", flush=True)

    tss.to_csv(tss_path, sep="\t", index=False, header=False)
    promoter_plus.to_csv(promoter_plus_path, sep="\t", index=False, header=False)
    promoter_minus.to_csv(promoter_minus_path, sep="\t", index=False, header=False)

    print(f"Plus-strand promoters: {len(promoter_plus)}", flush=True)
    print(f"Minus-strand promoters: {len(promoter_minus)}", flush=True)

    if overlaps_path.exists():
        print("Plotting Overlaps...", flush=True)
        plot_overlap_histograms(overlaps_path, histogram_dir)
        print(f"Saved overlap histograms to: {histogram_dir}", flush=True)
    else:
        print("promoter_overlaps_all.txt not found yet.", flush=True)


if __name__ == "__main__":
    main()