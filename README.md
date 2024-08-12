# KRANK shared scripts, results and dataset information

## Results and evaluation metrics
`results/cscores-10kSpecies_-_combined.csv` : taxonomic classification results for all tools on WoL queries
`results/all_tools-profiling_evaluation-CAMI1_hc.tsv` : abundance profiling metrics of all tools on CAMI-I high-complexity dataset
`results/cscores-10kSpecies_-_with_sizes.csv` : taxonomic classification results for all tools on WoL queries with library sizes of each tool
`results/resultsCAMI2-marine.tsv` : abundance profiling metrics of all tools on CAMI-II marine dataset
`results/cscores-10kSpecies_-_KRANK-candidates.csv` : comparison of different sizes and parameters for KRANK on WoL taxonomic classification
`results/cscores-10kSpecies_-_Kraken-II_4Gb.csv` : taxonomic classification results of Kraken-II on WoL using 4GB
`results/resultsCAMI2-strain_madness.tsv` : abundance profiling metrics of all tools on CAMI-II strain-madness dataset
`results/cscores-10kSpecies_-_CLARK.csv` : taxonomic classification results of CLARK on WoL
`results/cscores-10kSpecies_-_Kraken-II_16Gb.csv` : taxonomic classification results of Kraken-II on WoL using 16GB
`results/cscores-10kSpecies_-_CONSULT-II.csv` : taxonomic classification results of CONSULT on WoL using the default configuration
`results/running_times-query.tsv` : query and library construction running times
`results/cscores-10kSpecies_-_KRANK-rankingkmers_comparison.csv` : comparison of different heuristics for KRANK's selection on WoL taxonomic classification
`results/cscores-10kSpecies_-_Kraken-II_default.csv` : taxonomic classification results of Kraken-II on WoL using the default parameters
`results/cscores-10kSpecies_-_KRANK-sizeconst_comparison.csv` : comparison of different heuristics for KRANK's size constraint on WoL taxonomic classification

## Auxiliary data, dataset descriptions and taxonomy information
`data/ReferenceTaxonomy-nodes.dmp.gz` : WoL-v1 taxonomy nodes
`data/ref_taxa_counts.txt` : genome counts for each taxon in WoL-v1 with rank information
`data/10kBacteria-metadata.tsv` : WoL-v1 metadata including download links and additional information for reference genomes
`data/query_genomes_list.txt` : IDs of genomes used in read classification on WoL (download simulated reads [here](https://ter-trees.ucsd.edu/data/krank/KRANK-queries.tar.gz))
`data/ref_genome_counts` :  genome counts for each taxon in WoL-v1
`data/ReferenceTaxonomyRWoL-nodes.dmp.gz` : WoL-v1 taxonomy nodes reduced to species set
`data/taxonomy_lookup` : taxonomy lookup table used by CONSULT-II, parent list of each taxon
`data/dist_wrt_lastcommonrank.csv` : Jaccard similarity between randomly sampled genomes and their corresponding groups
`data/query_ranks.tsv` : taxonomy information for query genomes, ground truth for evaluation of taxonomic classification
`data/reference_genomes_list` : reference genomes and corresponding species
`data/uDance-ranks_tid.tsv` : WoL-v2 taxonomic ranks, some queries were retrieved from here
`data/10kBacteria-ranks_tid.tsv` : WoL-v1 taxonomic ranks, all genomes in the reference library
`data/dist_to_closest.txt` : closest reference genome of each query genome and their genomic distance similarity estimated by Mash
`data/download-links/all_download.txt` : all download links for WoL-v2 used in uDance
`data/download-links/download_final_extra_queries.txt` : download links for genomes that are not used in CONSULT-II paper
`data/download-links/genomes_uniq_uDance.txt`
`data/download-links/downloads-uDance_exc10k.txt`
`data/auxiliary/sampleg_dists.txt`
`data/auxiliary/dist-extra-to-closest.txt`
`data/auxiliary/uDance_exc10k-ranks_tid.tsv`
`data/auxiliary/uDance-genera_list.txt`
`data/auxiliary/uDance-species_list.txt`
`data/auxiliary/uDance_oneperfamily-ranks_tid.tsv`
`data/auxiliary/uDance_exc10k-order_info`
`data/auxiliary/closest_taxon_wrank.txt`
`data/auxiliary/dist-bacteria-to-closest.txt`
`data/auxiliary/uDance_exc10k-ranks_tid-downloadable.tsv`
`data/auxiliary/dist-to-closest.txt`
`data/auxiliary/dist-archaea-to-closest.txt`

## Scripts used to do empirical evaluation and create figures
`scripts/construct_taxonomy_lookup.py` : constructs the taxonomy lookup table for CONSULT-II from a taxonomy nodes file
`scripts/shrink_taxdump.py` : given taxonomy nodes and names files and a set of species, reduces taxonomy to the set of species of interest
`scripts/evaluate_CLARK.py` : custom script to evaluate the read classification output of CLARK, computes TP/FP/TN/FN for each rank and each read
`scripts/evaluate_KRANK.py` : custom script to evaluate the read classification output of KRANK, computes TP/FP/TN/FN for each rank and each read
`scripts/evaluate_CONSULTII.py` : custom script to evaluate the read classification output of CONSULT-II, computes TP/FP/TN/FN for each rank and each read
`scripts/evaluate_KrakenII.py` : custom script to evaluate the read classification output of Kraken-II, computes TP/FP/TN/FN for each rank and each read
`scripts/summarize_evaluations.py` : summarize TP/FN/TN/FN counts across ranks and genomes, should be used with the output of above evaluate_\*.py scripts
`scripts/prepprocess_methods_psummary.py` : uses distances in dist/dist_to_closest.txt to compute F1/precision/recall for different distance levels
`scripts/prepprocess_methods_summary.py` : uses distances in dist/dist_to_closest.txt to compute F1/precision/recall across different novelty bins
`scripts/prepprocess_methods_csummary.py` : uses distances in dist/dist_to_closest.txt to compute F1/precision/recall across taxon sizes
`scripts/match_closest_taxon.py`
`scripts/dist_wrt_lastcommonrank.py`
`scripts/get_taxa_count.py`
`scripts/count_taxa.sh`
`scripts/find_closest_taxon.sh`
`scripts/resource_benchmarking.R`
`scripts/shared_kmers_analysis.R`
`scripts/profiling_cami2_analysis.R`
`scripts/profiling_tool_comparision.R`
`scripts/comparison_-_withCONSULT-II.R`
`scripts/size_const_comparison-10kSpecies.R`
`scripts/kmer_ranking_comparison-10kSpecies.R`
`scripts/classification_comparison-10kSpecies.R`
`scripts/numgenomes_per_taxon-violinplot-10kSpecies.R`
`scripts/summary_analysis_cami2.R`
`scripts/weight_dist_simulations.R`
`scripts/query_info.R`

## Figures and illustrations
`figures/query_details.pdf`
`figures/profiling-cami2_combined-tool_comparison-l1_unifrac-avg_ranks.pdf`
`figures/classification_comparison-main-10kSpecies.pdf`
`figures/classification_comparison-wrt_memory.pdf`
`figures/classification_comparison-defaultsPrecisionRecall-10kSpecies.pdf`
`figures/classification_comparison-main_new-10kSpecies.pdf`
`figures/profiling-cami2_combined-tool_comparison-completeness_purity.pdf`
`figures/profiling-cami2_combined-tool_comparison-strain_madness.pdf`
`figures/size_const_comparison-wrt_group_size.pdf`
`figures/shared_kmers_portion.pdf`
`figures/profiling-cami2_combined-tool_comparison-l1_unifrac-avg_all.pdf`
`figures/classification_comparison-null_model-10kSpecies.pdf`
`figures/running_time-query.pdf`
`figures/classification_comparison-defaultsF1-10kSpecies.pdf`
`figures/classification_comparison_-_withCONSULT-II.pdf`
`figures/size_const_comparison-10kSpecies.pdf`
`figures/improvement_new_profiling.pdf`
`figures/expected_num_matches.pdf`
`figures/improvement_profiling-genome_size_correction.pdf`
`figures/shared_kmers_analysis.pdf`
`figures/num_genomes_per_taxon.pdf`
`figures/improvement_profiling-new_method.pdf`
`figures/classification_comparison-varying_memory-10kSpecies.pdf`
`figures/kmer_ranking_comparison-10kSpecies.pdf`
`figures/profiling-cami2_combined-tool_comparison-l1_unifrac.pdf`
`figures/krank-illustration.pdf`
`figures/krank-illustration.key`
`figures/profiling-cami2_combined-tool_comparison.pdf`
`figures/distance_to_closest.pdf`
`figures/profiling_tool_comparison.pdf`
