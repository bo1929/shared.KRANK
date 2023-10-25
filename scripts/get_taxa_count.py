import sys
import csv
from collections import defaultdict

if __name__ == "__main__":
    with open(sys.argv[1], 'r') as file:
        all_taxa_counts = map(lambda x: x.strip().split('\t'), file.readlines())
        all_taxa_counts = map(lambda x: (x[0], x[1:]), all_taxa_counts)
        all_taxa_counts = dict(map(lambda x: (x[0], list(map(lambda y: tuple(y.split(' ')), x[1]))), all_taxa_counts))
    with open(sys.argv[1], 'w') as file:
        file.write("Rank\tTaxon\tCount\n")
        for rank, tc_list in all_taxa_counts.items():
            for count, taxon in tc_list:
                file.write(f"{rank}\t{taxon}\t{count}\n")
