import sys
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from collections import defaultdict
import math

def parse_dist_file(filename, n_bins):
    sys.stderr.write("Start reading input file ...\n")
    processed_lines = 0
    occ = defaultdict(lambda: 0)
    maxval = 0
    for line in open(filename, 'r'):
        if processed_lines % 100000 == 0:
            sys.stderr.write("Processed " + str(processed_lines) + " lines.\n")
        fields = line.split()
        distance = int(fields[1])
        occ[distance] += 1
        if distance > maxval:
            maxval = distance
        processed_lines += 1

    sys.stderr.write("Create histogram ...\n")
    # produce a histogram with n_bin bins
    bin_width = max(1, math.ceil((maxval + 1) / n_bins))
    histogram = [0] * n_bins

    for dist, freq in occ.items():
        bin_idx = dist // bin_width
        histogram[bin_idx] += freq

    bin_edges = [int(i * bin_width) for i in range(n_bins + 1)]
    return histogram, bin_edges


def plot_histogram(histogram, bin_edges, filename):
    bin_width = bin_edges[1] - bin_edges[0] if len(bin_edges) > 1 else 1
    plt.figure(figsize=(10, 6))
    plt.bar(bin_edges[:len(histogram)], histogram,
            width=bin_width) #, align='edge', edgecolor='black')
    plt.xlabel("Distance")
    plt.ylabel("Count")
    plt.yscale('symlog')
    plt.title("Distance between unique syncmers")
    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close()


if __name__ == '__main__':

    if (len(sys.argv) < 3):
        sys.stderr.write("Usage: python3 plot_distance_stats.py <distances.tsv> <outfile.pdf> .\n")
        sys.exit(1)

    distance_file = sys.argv[1]
    outfile = sys.argv[2]
    n_bins = int(sys.argv[3])
    histogram, bin_edges = parse_dist_file(distance_file, n_bins)
    plot_histogram(histogram, bin_edges, outfile)
