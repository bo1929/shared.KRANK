paste -s -d'	' <(cut -f2 ../data/10kBacteria-ranks_tid.tsv | sort -n | uniq -c) \
  <(cut -f3 ../data/10kBacteria-ranks_tid.tsv | sort -n | uniq -c) \
  <(cut -f4 ../data/10kBacteria-ranks_tid.tsv | sort -n | uniq -c) \
  <(cut -f5 ../data/10kBacteria-ranks_tid.tsv | sort -n | uniq -c) \
  <(cut -f6 ../data/10kBacteria-ranks_tid.tsv | sort -n | uniq -c) \
  <(cut -f7 ../data/10kBacteria-ranks_tid.tsv | sort -n | uniq -c) \
  <(cut -f8 ../data/10kBacteria-ranks_tid.tsv | sort -n | uniq -c) \
  | tr -s ' ' | sed 's/	 /	/g' | sed 's/^ //g' | sed 's/^[^[:alpha:]]*//' > ../data/ref_taxa_counts.txt
python3 get_taxa_count.py ../data/ref_taxa_counts.txt
