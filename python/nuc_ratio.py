
import pandas as pd
import argparse
import scipy.stats as stats



def nucleosome_ratio(controls_bed, sample_bed):
    """Check the difference of cfDNA fragment length between sample and control data, by checking the mono-nucleosome and di-nucleosome fragment length 
    differences uses Wilcoxon rank sum test.
    
    """

    control_data = pd.read_csv(controls_bed, sep="\t", header=None)
    sample_data = pd.read_csv(sample_bed, sep="\t", header=None)


    # subtract the values in the third column from the values in the second column
    control_data['length'] = control_data.iloc[:, 2] - control_data.iloc[:, 1]

    sample_data['length'] = sample_data.iloc[:, 2] - sample_data.iloc[:, 1]


    mono_nuc_control = control_data[control_data['length'] < 150]
    mono_nuc_sample = sample_data[sample_data['length'] < 150]

    statistic_mono, pvalue_mono = stats.ranksums(mono_nuc_control, mono_nuc_sample)


    di_nuc_control = control_data[(control_data['length'] >= 275) & (control_data['length'] <= 400)]
    
    di_nuc_sample = sample_data[(sample_data['length'] >= 275) & (sample_data['length'] <= 400)]

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
