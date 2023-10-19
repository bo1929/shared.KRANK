import sys
import pandas as pd
from collections import defaultdict

r_method = pd.read_csv(f"../{sys.argv[3]}/{sys.argv[2]}-summary_evaluation-c{sys.argv[1]}.tsv", sep="\t", header=None)


def compute_metrics(df, name=""):
    scores = defaultdict(list)
    for bin_rank, sub_df in df.groupby([0, 1]):
        values_TP = sub_df.loc[sub_df[2] == "TP"][3].values
        values_TN = sub_df.loc[sub_df[2] == "TN"][3].values
        values_FP = sub_df.loc[sub_df[2] == "FP"][3].values
        values_FN = sub_df.loc[sub_df[2] == "FN"][3].values
        if values_TP.size > 0:
            count_TP = values_TP[0]
        else:
            count_TP = 0
        if values_TN.size > 0:
            count_TN = values_TN[0]
        else:
            count_TN = 0
        if values_FP.size > 0:
            count_FP = values_FP[0]
        else:
            count_FP = 0
        if values_FN.size > 0:
            count_FN = values_FN[0]
        else:
            count_FN = 0
        precision = (count_TP) / (count_TP+count_FP+10**-5)
        recall = (count_TP) / (count_TP+count_FN+10**-5)
        scores["Taxonomic_Rank"].append(bin_rank[1])
        scores["Distance_to_closest"].append(f"[{eval(bin_rank[0])[0]}, {eval(bin_rank[0])[1]})")
        scores["Precision"].append(precision)
        scores["Recall"].append(recall)
        scores["F1"].append((2*precision*recall)/(precision+recall+10**-5))
        scores["Method"].append(name)
    return pd.DataFrame(scores)

scores_method = compute_metrics(r_method, f"{sys.argv[2]} (0.{sys.argv[1]})")
scores_method.to_csv(f"../{sys.argv[3]}/{sys.argv[2]}-summary_scores_methods-{sys.argv[1]}.csv")
