import sys
import collections as coll

BINS = [(0.0, 0.001), (0.001, 0.02), (0.02, 0.06), (0.06, 0.12), (0.12, 0.16), (0.16,0.22), (0.22, 0.35)]
CONT = bool(int(sys.argv[3]))


def get_bin(dist):
    bin_i = (0.35, 1)
    for bin_i in BINS:
        if (dist >= bin_i[0]) and (dist < bin_i[1]):
            break
    return bin_i

if __name__ == "__main__":
    results = coll.defaultdict(lambda: coll.defaultdict(lambda: coll.defaultdict(int)))

    with open(sys.argv[2], "r") as dist_file:
        dist_map = {genome[:-4]:dist for genome, dist in map(lambda x: (x[0], x[2]), map(lambda x: x.strip().split("\t"), dist_file.readlines()))}

    if CONT:
        with open(sys.argv[1], "r") as results_file:
            for line in results_file.readlines()[1: ]:
                ls = line.strip().split(",")
                results[ls[1]]["species"][ls[5]] +=1
                results[ls[1]]["genus"][ls[6]] +=1
                results[ls[1]]["family"][ls[7]] +=1
                results[ls[1]]["order"][ls[8]] +=1
                results[ls[1]]["class"][ls[9]] +=1
                results[ls[1]]["phylum"][ls[10]] +=1
                results[ls[1]]["superkingdom"][ls[11]] +=1

        for genome, results_bin in results.items():
            for key_rank, results_rank in results_bin.items():
                for key_type, count in results_rank.items():
                    dist_to_closest = float(dist_map[genome])
                    print(f"{genome}\t{dist_to_closest}\t{key_rank}\t{key_type}\t{count}")
    else:
        with open(sys.argv[1], "r") as results_file:
            for line in results_file.readlines()[1: ]:
                ls = line.strip().split(",")
                dist_to_closest = float(dist_map[ls[1]])
                bin_genome = get_bin(dist_to_closest)
                results[bin_genome]["species"][ls[5]] +=1
                results[bin_genome]["genus"][ls[6]] +=1
                results[bin_genome]["family"][ls[7]] +=1
                results[bin_genome]["order"][ls[8]] +=1
                results[bin_genome]["class"][ls[9]] +=1
                results[bin_genome]["phylum"][ls[10]] +=1
                results[bin_genome]["superkingdom"][ls[11]] +=1

        for key_bin, results_bin in results.items():
            for key_rank, results_rank in results_bin.items():
                for key_type, count in results_rank.items():
                    print(f"{key_bin}\t{key_rank}\t{key_type}\t{count}")
