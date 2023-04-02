
import pandas as pd
import argparse
import scipy.stats as stats



def function(controls_bed):
    """
    check the cfDNA data that have lengths above 200 bases
    """
    control_data = pd.read_csv(controls_bed, sep="\t", header=None)


    # subtract the values in the third column from the values in the second column
    control_data['length'] = control_data.iloc[:, 2] - control_data.iloc[:, 1]



    di_nuc_control = control_data[(control_data['length'] >= 200)]
    print(di_nuc_control.head())




def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('input', help='control file')

    args = parser.parse_args()

    function(args.input)
    

if __name__ == '__main__':
    main()
