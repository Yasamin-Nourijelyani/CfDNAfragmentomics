# Write some python boiler plate code
# Import the necessary modules
from tqdm import tqdm
import numpy as np
import argparse
import csv
from collections import defaultdict
import matplotlib.pyplot as plt
from scipy import signal

def plot_coverage(coverage, name):
    """
    Plot the coverage of the nucleosome occupancy in a line plot
    """
    x = np.linspace(-1000, 1000, 2001)
    # smooth with a savitzky golay filter
    y = signal.savgol_filter(coverage, 51, 3)
    plt.plot(x, y)
    plt.savefig(f'{name}.NO.png')


def aggregate_coverage(coverage, tfbs, window):
    """
    aggregate the coverage vectors across every window in the tfbs file

    :param coverage: a dictonary of chromosome names and coverage vectors
    :paramtfbs: path to the tfbs file
    :param window: the size of the window to aggregate over
    :return: aggregated coverage vector
    """
    aggregate_coverage = np.zeros(2*window+1)
    with open(tfbs, 'r') as tfbs:
        reader = csv.reader(tfbs, delimiter='\t')
        header = next(reader)
        # map the header to the index of the column
        header = {header[i]: i for i in range(len(header))}
        print("Aggregating coverage")
        for row in tqdm(reader):
            chr = row[header['chr']]
            midpoint = int(row[header['mid_position']])
            start = midpoint - window
            end = midpoint + window
            aggregate_coverage += coverage[chr][start:end+1]

    return aggregate_coverage

def get_coverage(bed_file):
    """
    Get the coverage at every locus in the interval [start, end]

    :param bed_file: The bed file to get the coverage from
    :return: A list of coverage values
    """

    coverage = defaultdict(lambda: np.zeros(249*10**6))
    with open(bed_file, 'r') as bed:
        reader = csv.reader(bed, delimiter='\t')
        print("Getting coverage")
        for row in tqdm(reader):
            chr = row[0]
            start = int(row[1])
            end = int(row[2])
            for i in range(start, end+1):
                coverage[chr][i] += 1

    return coverage



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input bed file")
    parser.add_argument('-t', '--tfbs', help='TFBS file')
    parser.add_argument('-w', '--window', help='Window size to search around TFBS', type=int)
    args = parser.parse_args()

    coverage = get_coverage(args.input)
    aggregated_coverage = aggregate_coverage(coverage, args.tfbs, args.window)
    plot_coverage(aggregated_coverage, args.input.split('.')[1])

if __name__ == '__main__':
    main()
