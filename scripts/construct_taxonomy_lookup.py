import csv
import argparse
from pathlib import Path
from collections import defaultdict

parser = argparse.ArgumentParser(description="Creates a taxonomic ancestor lookup table.")
parser.add_argument(
    "--input-path",
    type=str,
    required=True,
)
parser.add_argument(
    "--output-dir",
    type=str,
    required=False,
    default="./taxonomy_lookup/"
)
args = parser.parse_args()


PATH_RANKS = Path(args.input_path)
DIR_LOOKUPS = Path(args.output_dir)
DIR_LOOKUPS.mkdir(parents=True, exist_ok=True)

taxa_levels = ["species", "genus", "family", "order", "class", "phylum", "kingdom"]
taxa_level_ranks = {
    "kingdom": 7,
    "phylum": 6,
    "class": 5,
    "order": 4,
    "family": 3,
    "genus": 2,
    "species": 1,
}

if __name__ == "__main__":
    with open(PATH_RANKS, 'r', newline='') as csvfile:
        reader = list(csv.DictReader(csvfile, delimiter="\t"))
    ancestor_lookup = defaultdict(list)
    genome_lookup = defaultdict(int)
    level_lookup = defaultdict(str)

    for row_ranks in reader:
        genome_lookup[row_ranks["genome"]] = row_ranks["species"]
        for ix, taxa1 in enumerate(taxa_levels[:]):
            tid1 = row_ranks[taxa1]
            if tid1 not in ancestor_lookup.keys() and tid1 != 0:
                for taxa2 in taxa_levels[ix:]:
                    ancestor_lookup[tid1].append(row_ranks[taxa2])
                level_lookup[row_ranks[taxa1]] = taxa1
            else:
                break

    with open(DIR_LOOKUPS / "ancestor_lookup", "w") as f:
        lookup_str = ""
        for taxa, ancestors in ancestor_lookup.items():
            lookup_str += f"{taxa} {','.join([str(a_taxa) for a_taxa in ancestors])}"
            lookup_str += "\n"
        f.write(lookup_str)

    with open(DIR_LOOKUPS / "genome_lookup", "w") as f:
        lookup_str = ""
        for key, val in genome_lookup.items():
            lookup_str += f"{key} {val}"
            lookup_str += "\n"
        f.write(lookup_str)

    with open(DIR_LOOKUPS / "level_lookup", "w") as f:
        lookup_str = ""
        for key, val in level_lookup.items():
            lookup_str += f"{key} {val}"
            lookup_str += "\n"
        f.write(lookup_str)
