#!/usr/bin/env python3
# usage: plot_nx.py distances.tsv nx_plot.png
import sys
from itertools import groupby, accumulate
from bisect import bisect_left
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

if len(sys.argv) < 3:
    sys.stderr.write("Error: too few arguments.\n")
    sys.stderr.write("Usage: python3 plot_nx_curve.py <input> <output>\n")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file) as f:
    f.readline()
    rows = (line.split("\t") for line in f)
    for src, grp in groupby(rows, key=lambda t: t[0]):
        pieces = [(int(t[1]), int(t[2])) for t in grp]
        d = [p[0] for p in pieces]
        cum = list(accumulate(p[0] * p[1] for p in pieces))
        total = cum[-1]
        xs = [i / 10 for i in range(1, 1001)]
        ys = [d[min(bisect_left(cum, x / 100 * total), len(d) - 1)] for x in xs]
        plt.step(xs, ys, where="post", color="C0", alpha=0.15, lw=0.6)

plt.yscale("log")
plt.xlabel("x (%)")
plt.title("Nx")
plt.ylabel("distances between unique syncmers (bp)")
plt.tight_layout()
plt.savefig(output_file, dpi=200)
