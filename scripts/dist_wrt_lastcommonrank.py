import pandas as pd
import numpy as np
from collections import defaultdict

ranks = pd.read_csv("../data/10kBacteria-ranks_tid.tsv", sep='\t')
ranks.index = ranks["genome"]
with open("../data/auxiliary/sampleg_dists.txt", 'r') as f:
    all_dists = list(map(lambda x: x.strip().split(), f.readlines()))

dd = defaultdict(list)
for pair in all_dists:
    g1, g2 = pair[0], pair[1]
    d = float(pair[2])
    eq = np.logical_and((ranks.loc[g1] == ranks.loc[g2]), ranks.loc[g2])
    ll = np.nonzero(eq)[0]
    if len(ll) > 0:
        match_rank = (list(ranks.columns)[max(ll)])
    else:
        match_rank = "root"
    if d < 1:
        dd["g1"].append(g1)
        dd["g2"].append(g1)
        dd["d"].append(d)
        dd["rank"].append(match_rank)
        dd["JI"].append(eval(pair[-1]))

xd = pd.DataFrame(dd)

xd.sample(n=500000).to_csv("../data/dist_wrt_lastcommonrank.csv")
