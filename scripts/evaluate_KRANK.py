import argparse
import multiprocessing
import pandas as pd

from pathlib import Path
from collections import defaultdict

parser = argparse.ArgumentParser(description="Evaulates performance of CONSULT-II.")
parser.add_argument(
    "--results-path",
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
parser.add_argument(
    "--tvth",
    type=float,
    required=False,
    default=0.00
)
args = parser.parse_args()

results_path = Path(args.results_path)
QUERY_RANKS_PATH =  Path(args.query_ranks_path)
REFERENCE_RANKS_PATH = Path(args.reference_ranks_path)
TAXONOMY_DATABASE_PATH = Path(args.taxonomy_database_path)
OUTPUT_DIR = Path(args.output_dir)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
tvth = args.tvth

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

    # example filename : output_1000x-G000251165.txt
    genome = result_file.stem.split("-")[-1]
    true_taxa = query_ranks.loc[genome]

    with open(result_file) as f:
        all_reads = f.readlines()
    for read_line in all_reads:
        read_line_cols = read_line.split("\t")
        read_name = read_line_cols[0]
        read_form = read_line_cols[1]
        pred_tID = int(read_line_cols[2])
        pred_vote = float(read_line_cols[3])
        max_vote = float(read_line_cols[4])

        evaluation_dict["genome"].append(genome)
        evaluation_dict["read"].append(read_name)

        pred_taxa = {}
        rank = None

        if pred_vote > tvth:
            evaluation_dict["vote"].append(pred_vote)
            evaluation_dict["voteNormalized"].append(pred_vote / (max_vote + 10**(-6)))

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
        else:
            evaluation_dict["vote"].append(pred_vote)
            evaluation_dict["voteNormalized"].append(pred_vote / (max_vote + 10**(-6)))

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
    # Clean whitespaces and newlines etc.
    reference_taxonomy = reference_taxonomy.applymap(lambda x: x.strip() if isinstance(x, str) else x)

    # Make sure we use superkingdom consistently.
    reference_ranks = reference_ranks.rename(columns={"kingdom": "superkingdom"})
    query_ranks = query_ranks.rename(columns={"kingdom": "superkingdom"})

    # Set DataFrame indices to genome names.
    query_ranks.index = query_ranks["genome"]
    reference_ranks.index = reference_ranks["genome"]
    reference_ranks = reference_ranks.drop("genome", axis=1)
    query_ranks = query_ranks.drop("genome", axis=1)

    # Set DataFrame indices to taxonomic IDs.
    reference_taxonomy.index = reference_taxonomy[0]
    df = evaluate_classification(results_path)
    genome = results_path.stem.split("-")[-1]
    str_tvth = "{:.2f}".format(tvth)
    df.to_csv(OUTPUT_DIR / f"{genome}-evaluations-c{str_tvth[2:]}.csv")
