"""
07_plot_hsf1_binding_distribution.py

Plots the genomic distribution of HSF-1 binding sites from FIMO output.

Input:  results/fimo_out_dedupes/fimo.tsv
Output: results/plots/hsf1_distribution/
    - chrI.pdf ... chrX.pdf  (per-chromosome histograms)
    - hsf1_all_overlaid.pdf  (all chromosomes overlaid)
    - hsf1_sorted_by_chr_position.csv
"""

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import seaborn as sns

CHROM_ORDER = ["N2_chrI", "N2_chrII", "N2_chrIII", "N2_chrIV", "N2_chrV", "N2_chrX"]
CHROM_COLORS = ["#4ec9b0", "#569cd6", "#c586c0", "#ce9178", "#dcdcaa", "#f44747"]


def load_fimo(input_path: Path) -> pd.DataFrame:
    df = pd.read_csv(input_path, sep="\t", comment="#")
    df["chromosome"] = df["sequence_name"].str.split("|").str[-1].str.extract(
        r"(N2_chr[IVX]+)"
    )[0]
    df["chromosome"] = pd.Categorical(df["chromosome"], categories=CHROM_ORDER, ordered=True)
    n_before = len(df)
    df = df.drop_duplicates()
    n_dropped = n_before - len(df)
    if n_dropped:
        print(f"Dropped {n_dropped} duplicate rows ({n_before} → {len(df)})")
    return df.sort_values(["chromosome", "start"]).reset_index(drop=True)


def apply_style() -> None:
    sns.set_style("dark")
    plt.rcParams.update({
        "figure.facecolor":  "black",
        "axes.facecolor":    "#111111",
        "text.color":        "white",
        "axes.labelcolor":   "white",
        "xtick.color":       "white",
        "ytick.color":       "white",
        "axes.spines.top":   False,
        "axes.spines.right": False,
    })


def format_axis(ax: plt.Axes) -> None:
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x/1e7:.1f}"))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(2e6))


def plot_individual(df: pd.DataFrame, outdir: Path) -> None:
    for chrom in CHROM_ORDER:
        group = df[df["chromosome"] == chrom]
        if group.empty:
            print(f"Warning: no data for {chrom}, skipping.")
            continue

        fig, ax = plt.subplots(figsize=(8, 4))
        fig.patch.set_facecolor("black")
        ax.hist(group["start"], bins=200, color="green", alpha=0.5, linewidth=0)

        chrom_label = chrom.replace("N2_", "")
        ax.set_title(chrom_label, color="white", fontsize=14)
        ax.set_xlabel("Genomic Position (x10^7)", color="white", fontsize=11)
        ax.set_ylabel("Frequency", color="white", fontsize=11)
        format_axis(ax)
        plt.tight_layout()

        out = outdir / f"{chrom_label}.pdf"
        plt.savefig(out, bbox_inches="tight", facecolor="black")
        plt.close()
        print(f"Saved {out}")


def plot_overlaid(df: pd.DataFrame, outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 5))
    fig.patch.set_facecolor("black")

    for chrom, color in zip(CHROM_ORDER, CHROM_COLORS):
        group = df[df["chromosome"] == chrom]
        ax.hist(
        group["start"],
        bins=200,
        color=color,
        alpha=0.5,
        label=chrom.replace("N2_", ""),
        linewidth=0,
    )

    ax.set_title("HSF-1 binding sites — all chromosomes", color="white", fontsize=14)
    ax.set_xlabel("Genomic Position (x10^7)", color="white", fontsize=11)
    ax.set_ylabel("Frequency", color="white", fontsize=11)
    format_axis(ax)
    ax.legend(frameon=False, labelcolor="white", fontsize=9)
    plt.tight_layout()

    out = outdir / "hsf1_all_overlaid.pdf"
    plt.savefig(out, bbox_inches="tight", facecolor="black")
    plt.close()
    print(f"Saved {out}")


def main() -> None:
    root = Path(__file__).resolve().parents[1]

    input_path = root / "results" / "fimo_out_dedupes" / "fimo.tsv"
    outdir = root / "results" / "plots" / "hsf1_distribution"
    csv_out = root / "data" / "processed" / "hsf1_sorted_by_chr_position.csv"

    outdir.mkdir(parents=True, exist_ok=True)

    print("Loading FIMO results...", flush=True)
    df = load_fimo(input_path)
    print(f"Loaded {len(df)} binding sites across {df['chromosome'].nunique()} chromosomes.")

    df.to_csv(csv_out, index=False)
    print(f"Saved sorted CSV: {csv_out}")

    apply_style()
    plot_individual(df, outdir)
    plot_overlaid(df, outdir)


if __name__ == "__main__":
    main()
