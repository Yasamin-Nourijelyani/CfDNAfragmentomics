
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
    
    Sample use:
      
      python ../python/nuc_ratio.py --sample lung_cancer_data/EE88183.hg38.frag.bed healthy_data/EE86217.hg38.frag.bed
    
    """
    # all of the length data
    lengths_control = []

    with open(controls_bed, 'r') as c_bed:
        reader = csv.reader(c_bed, delimiter='\t')
        print("Getting sample lengths")
        for row in tqdm(reader):
            start = int(row[1])
            end = int(row[2])
            length = end - start
            lengths_control.append(length)

    
    lengths_sample = []
    with open(sample_bed, 'r') as s_bed:
        reader = csv.reader(s_bed, delimiter='\t')
        print("Getting sample lengths")
        for row in tqdm(reader):
            start = int(row[1])
            end = int(row[2])
            length = end - start
            lengths_sample.append(length)



    mono_nuc_control = lengths_control[lengths_control < 150]
    mono_nuc_sample = lengths_sample[lengths_sample < 150]

    statistic_mono, pvalue_mono = stats.ranksums(mono_nuc_control, mono_nuc_sample)


    di_nuc_control = lengths_control[(lengths_control >= 275) & (lengths_control <= 400)]
    
    di_nuc_sample = lengths_sample[(lengths_sample >= 275) & (lengths_sample <= 400)]

    statistic_di, pvalue_di = stats.ranksums(di_nuc_control, di_nuc_sample)




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
