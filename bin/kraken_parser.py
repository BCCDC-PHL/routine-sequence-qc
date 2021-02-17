#!/usr/bin/env python3

import argparse
import json
import sys


def main(args):

    taxonomic_levels = {'U', 'D', 'K', 'P', 'C', 'O', 'F', 'G', 'S'}
    # assert args.taxonomic_level in rank_codes, "Rank must be one of [(U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies.]"

    total_reads = 0
    unclassified_reads = 0
    total_reads_reported = 0
    other_reads = 0

    headers = [
        "percent_reads_in_clade",
        "num_reads_in_clade",
        "ncbi_taxonomy_id",
        "taxon_name",
    ]

    print('\t'.join(headers))

    with open(args.kraken_report, 'r') as f:
        for line in f:
            record = {}
            split_line = [record.strip() for record in line.strip().split("\t")]

            record['percent_reads_clade'] = float(split_line[0])
            record['num_reads_clade'] = int(split_line[1])
            record['num_reads_taxon'] = int(split_line[2])
            record['taxonomic_level'] = split_line[3]
            record['ncbi_taxonomy_id'] = split_line[4]
            record['clade_name'] = split_line[5]

            if record['clade_name'] == "unclassified":
                unclassified_reads = record['num_reads_clade']
                total_reads += record['num_reads_clade']
            elif record['clade_name'] == "root":
                total_reads += record['num_reads_clade']

            # record['percent_reads_in_clade'] = record['num_reads_clade'] / float(total_reads) * 100

            if record['taxonomic_level'] == args.taxonomic_level and record['percent_reads_clade'] >= args.threshold_percent:
                print('\t'.join([str(record['percent_reads_clade']), str(record['num_reads_clade']), str(record['ncbi_taxonomy_id']), record['clade_name']]))
                total_reads_reported += record['num_reads_clade']


    print(str('%.3f' % (unclassified_reads / float(total_reads) * 100)) + '\t' + str(unclassified_reads) + '\t\t' + 'unclassified')
    total_reads_reported += unclassified_reads
    other_reads = total_reads - total_reads_reported
    print(str('%.3f' % (other_reads / float(total_reads) * 100)) + '\t' + str(other_reads) + '\t\t' + 'other')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('kraken_report')
    parser.add_argument('-l', '--taxonomic_level')
    parser.add_argument('-t', '--threshold_percent', type=float)
    args = parser.parse_args()
    main(args)
