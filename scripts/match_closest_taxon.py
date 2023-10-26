import sys

if __name__ == "__main__":
    with open(sys.argv[1], 'r') as inf:
        r1 = map(lambda x: x.strip().split(), inf.readlines())
    with open(sys.argv[2], 'r') as inf:
        r2 = map(lambda x: x.strip().split(), inf.readlines())
    for pair in zip(r1, r2):
        last = 1
        c = 0
        for i in zip(*pair):
            c+=1
            if c > 1:
                if i[0] != i[1]:
                    break
                else:
                    if i[0] != '0':
                         last = i[0]
        print(last)
