
import pandas as pd
import argparse
import scipy.stats as stats
from tqdm import tqdm
import numpy as np
import argparse
import csv
from collections import defaultdict
import matplotlib.pyplot as plt
from scipy import signal



def nucleosome_ratio(controls_bed, sample_bed):
    """Check the difference of cfDNA fragment length between sample and control data, by checking the mono-nucleosome and di-nucleosome fragment length 
    differences uses Wilcoxon rank sum test.
    
    """
    # mononucleosome length data
    lengths_control_mono = []

    lengths_control_di = []


    with open(controls_bed, 'r') as c_bed:
        reader = csv.reader(c_bed, delimiter='\t')
        print("Getting control lengths")
        for row in tqdm(reader):
            start = int(row[1])
            end = int(row[2])
            length = end - start
            if length < 150:
                lengths_control_mono.append(length)
            elif length >= 275 & length <= 400:
                lengths_control_di.append(length) 

    
    lengths_sample_mono = []
    lengths_sample_di = []

    with open(sample_bed, 'r') as s_bed:
        reader = csv.reader(s_bed, delimiter='\t')
        print("Getting sample lengths")
        for row in tqdm(reader):
            start = int(row[1])
            end = int(row[2])
            length = end - start
            if length < 150:
                lengths_sample_mono.append(length)
            elif length >= 275 & length <= 400:
                lengths_sample_di.append(length)




    statistic_mono, pvalue_mono = stats.ranksums(lengths_control_mono, lengths_sample_mono)



    statistic_di, pvalue_di = stats.ranksums(lengths_control_di, lengths_sample_di)

    # the mean of the control samples for mono-num
    mean_mono_c = sum(lengths_control_mono)/len(lengths_control_mono)

     # the mean of the control samples for mono-num
    mean_di_c = sum(lengths_control_di)/len(lengths_control_di)

        # the mean of the control samples for mono-num
    mean_mono_s = sum(lengths_sample_mono)/len(lengths_sample_mono)

        # the mean of the control samples for mono-num
    mean_di_s = sum(lengths_sample_di)/len(lengths_sample_di)

    # are the lengths of control mono-nuc greater than sample mono-nuc?
    if mean_mono_c >= mean_mono_s:
        print("length of mono-nuc is greater in the control compared to sample - if significantly differnt, this indicates cancer")
    else:
        print("length of mono-nuc is smaller in the control compared to sample")

    if mean_di_c >= mean_di_s:
        print("length of di-nuc is greater in the control compared to sample - if significantly differnt, this indicates cancer")
    else:
        print("length of mono-nuc is smaller in the control compared to sample")



    # print the resulting DataFrame
    print(f"The mononucleosome Wilcox statistic: {statistic_mono}")

    print(f"The mononucleosome Wilcox p-value: {pvalue_mono}")

    print(f"The dinucleosome Wilcox statistic: {statistic_di}")

    print(f"The dinucleosome Wilcox p-value:{pvalue_di}")


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('input', help='control file')
    parser.add_argument('-s', '--sample', help='sample file')

    args = parser.parse_args()

    nucleosome_ratio(args.input, args.sample)
    

if __name__ == '__main__':
    main()
