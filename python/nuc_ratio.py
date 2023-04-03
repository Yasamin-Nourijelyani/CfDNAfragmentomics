
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




    # print the resulting DataFrame
    print(statistic_mono)

    print(pvalue_mono)

    print(statistic_di)

    print(pvalue_di)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('input', help='control file')
    parser.add_argument('-s', '--sample', help='sample file')

    args = parser.parse_args()

    nucleosome_ratio(args.input, args.sample)
    

if __name__ == '__main__':
    main()
