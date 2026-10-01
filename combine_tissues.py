#!/usr/bin/env python

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Concatenate TSV files with tissue labels.')
    parser.add_argument('tsv_files', nargs='+', help='Input TSV files.')
    parser.add_argument('-tissues', nargs='+', required=True, help='List of tissue labels, same length as tsv_files.')
    parser.add_argument('-out', required=True, help='Output filename.')

    args = parser.parse_args()

    if len(args.tsv_files) != len(args.tissues):
        parser.error('The number of tissue labels must match the number of input TSV files.')

    dataframes = []
    print("Summary:")
    for tsv_file, tissue_label in zip(args.tsv_files, args.tissues):
        df = pd.read_csv(tsv_file, sep='\t')
        df['tissue'] = tissue_label
        num_rows = len(df)
        print(f"File: {tsv_file}, Tissue: {tissue_label}, Rows: {num_rows}")
        dataframes.append(df)

    concatenated_df = pd.concat(dataframes, ignore_index=True)
    concatenated_df.to_csv(args.out, sep='\t', index=False)

if __name__ == '__main__':
    main()
