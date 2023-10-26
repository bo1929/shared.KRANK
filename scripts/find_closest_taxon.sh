cut -f1 ../data/dist_to_closest.txt | xargs -I{} grep -w {} ../data/uDance-ranks_tid.tsv > q1
cut -f2 ../data/dist_to_closest.txt | xargs -I{} grep -w {} ../data/uDance-ranks_tid.tsv > q2
paste <(python match_closest_taxon.py q1 q2) <(cut -f3 ../data/dist_to_closest.txt) > ../data/closest_taxon.txt && rm q1 q2
