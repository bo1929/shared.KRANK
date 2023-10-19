import argparse
import multiprocessing
import pandas as pd

from pathlib import Path
from collections import defaultdict

parser = argparse.ArgumentParser(description="Evaulates performance of Kraken2.")
parser.add_argument(
    "--results-dir",
    type=str,
    required=True,
)
parser.add_argument(
    "--query-ranks-path",
    type=str,
    required=True,
)
parser.add_argument(
    "--reference-ranks-path",
    type=str,
    required=True,
)
parser.add_argument(
    "--taxonomy-database-path",
    type=str,
    required=True,
)
parser.add_argument(
    "--output-dir",
    type=str,
    required=False,
    default="./"
)
args = parser.parse_args()

NUM_THREADS = 12

RESULTS_DIR = Path(args.results_dir)
QUERY_RANKS_PATH =  Path(args.query_ranks_path)
REFERENCE_RANKS_PATH = Path(args.reference_ranks_path)
TAXONOMY_DATABASE_PATH = Path(args.taxonomy_database_path)
OUTPUT_DIR = Path(args.output_dir)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

taxa_order = list(
    ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]
)
taxa_order_rank = {
    "superkingdom": 6,
    "phylum": 5,
    "class": 4,
    "order": 3,
    "family": 2,
    "genus": 1,
    "species": 0,
}

def evaluate_classification(result_file):
    evaluation_dict = defaultdict(list)

    # example filename : output_1000x_G000251165.txt
    genome = result_file.stem.split("_")[-1]
    true_taxa = query_ranks.loc[genome]

    with open(result_file) as f:
        all_reads = f.readlines()
    for read_line in all_reads:
        read_line_cols = read_line.split("\t")
        read_name = read_line_cols[1]
        is_classified = read_line_cols[0]
        pred_tID = int(read_line_cols[2])

        pred_taxa = {}
        rank = None

        if is_classified != "U":
            parent, rank = (
                reference_taxonomy.at[pred_tID, 1],
                reference_taxonomy.at[pred_tID, 2],
            )
            pred_taxa[rank] = pred_tID

            while parent != 1:
                pred_tID = parent
                parent, rank = (
                    reference_taxonomy.at[parent, 1],
                    reference_taxonomy.at[parent, 2],
                )
                pred_taxa[rank] = pred_tID

        evaluation_dict["genome"].append(genome)
        evaluation_dict["read"].append(read_name)

        seen = False
        miss = False

        for idx, rank in enumerate(taxa_order):
            if true_taxa[rank] == 0:
                if rank != "species":
                    evaluation_dict[rank].append(evaluation_dict[taxa_order[idx - 1]][-1])
                else:
                    raise ValueError("Taxonomic ID is 0 at species level.")
            elif rank in pred_taxa.keys():
                if pred_taxa[rank] == true_taxa[rank]:
                    for i in range(idx, len(taxa_order)):
                        evaluation_dict[taxa_order[i]].append("TP")
                    break
                else:
                    evaluation_dict[rank].append("FP")
                seen = True
            elif not seen:
                if (not miss) and (true_taxa[rank] != reference_ranks[rank]).all():
                    evaluation_dict[rank].append("TN")
                else:
                    evaluation_dict[rank].append("FN")
                    miss = True
            else:
                if rank != "species":
                    evaluation_dict[rank].append(evaluation_dict[taxa_order[idx - 1]][-1])
                else:
                    raise ValueError("Some unkown condition or conflicting taxonomic ID.")
    return pd.DataFrame(evaluation_dict)


if __name__ == "__main__":
    query_ranks = pd.read_csv(QUERY_RANKS_PATH, sep="\t")
    reference_ranks = pd.read_csv(REFERENCE_RANKS_PATH, sep="\t")
    reference_taxonomy = pd.read_csv(
        TAXONOMY_DATABASE_PATH,
        sep="|",
        header=None,
        skipinitialspace=True
    )
    reference_taxonomy = reference_taxonomy.applymap(lambda x: x.strip() if isinstance(x, str) else x)

    kraken_result_dir = Path(RESULTS_DIR)
    # Note that we have Kraken2 output files with filenames: "output_1000x_<GENOME_NAME>.txt".

    # Make sure we use superkingdom consistently.
    reference_ranks = reference_ranks.rename(columns={"kingdom": "superkingdom"})
    query_ranks = query_ranks.rename(columns={"kingdom": "superkingdom"})

    # Clean whitespaces and newlines etc.

    # Set DataFrame indices to genome names.
    query_ranks.index = query_ranks["genome"]
    reference_ranks.index = reference_ranks["genome"]
    reference_ranks = reference_ranks.drop("genome", axis=1)
    query_ranks = query_ranks.drop("genome", axis=1)

    # Set DataFrame indices to taxonomic IDs.
    reference_taxonomy.index = reference_taxonomy[0]

    pool = multiprocessing.Pool(NUM_THREADS)
    all_evaluations = [
        pool.map(
            evaluate_classification,
            [
                file
                for file in list(kraken_result_dir.iterdir())
                if not file.stem.startswith(".")
            ],
        )
    ]
    pd.concat(all_evaluations[0]).to_csv(OUTPUT_DIR / "Kraken2-eval.csv")
