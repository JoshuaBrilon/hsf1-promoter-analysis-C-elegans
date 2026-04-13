from pathlib import Path
import pandas as pd


def main():
    root = Path(__file__).resolve().parents[1]

    input_path = root / "results" / "fimo_out_dedupes" / "fimo.tsv"
    output_path = root / "data" / "processed" / "fimo_gene_summary.txt"

    print("Loading FIMO results...", flush=True)
    fimo = pd.read_csv(input_path, sep="\t", comment="#")

    # Expecting sequence_name headers like:
    # gene_name|N2_chrI:10497-11497(+)
    fimo["gene_name"] = fimo["sequence_name"].str.split("|", n=1).str[0]

    summary = (
        fimo.groupby("gene_name")
        .agg(
            num_fimo_hits=("sequence_name", "size"),
            num_unique_promoters_hit=("sequence_name", "nunique"),
            best_score=("score", "max"),
            best_p_value=("p-value", "min"),
            best_q_value=("q-value", "min"),
        )
        .reset_index()
    )

    summary = summary.sort_values(
        by=["num_fimo_hits", "best_q_value", "best_p_value"],
        ascending=[False, True, True]
    )

    summary.to_csv(output_path, sep="\t", index=False, float_format="%0.2e")

    print(f"Genes summarized: {len(summary)}", flush=True)
    print(f"Saved: {output_path}", flush=True)


if __name__ == "__main__":
    main()