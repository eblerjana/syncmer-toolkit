import sys
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

def parse_hist_file(filename):
    histogram = []
    line_idx = 0
    for line in open(filename, 'r'):
        fields = line.strip().split()
        assert int(fields[0]) == line_idx
        count = int(fields[1])
        line_idx += 1
        histogram.append(count)
    return histogram

def plot_histogram(histogram, filename):
    print('Plotting')
    with PdfPages(filename) as pdf:
        # first plot histogram
        labels = [str(i) for i in range(len(histogram))]
        counts = histogram
        x_values = [i for i in range(len(histogram))]

        fig, ax = plt.subplots()
        ax.bar(x_values, counts)
        plt.title("unique syncmers in K haplotypes")
        plt.yscale('symlog')
        ax.set_xlabel("K")
        ax.set_ylabel("Count")
        pdf.savefig()
        plt.close()

        # plot summed up histogram
        cum_hist = []
        sum = 0
        for h in histogram:
            sum += h
            cum_hist.append(sum)

        fig, ax = plt.subplots()
        ax.bar(x_values, cum_hist)
        plt.title("unique syncmers in <= K haplotypes")
        plt.yscale('symlog')
        ax.set_xlabel("K")
        ax.set_ylabel("Count")
        pdf.savefig()
        plt.close()




if __name__ == '__main__':

    if (len(sys.argv) < 3):
        sys.stderr.write("Usage: python3 plot_syncmer_stats.py <histogram.tsv> <outfile.pdf> .\n")
        sys.exit(1)

    histogram_file = sys.argv[1]
    outfile = sys.argv[2]
    histogram = parse_hist_file(histogram_file)
    plot_histogram(histogram, outfile)
