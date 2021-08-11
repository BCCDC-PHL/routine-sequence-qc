#!/usr/bin/env python

import argparse
import csv
import json

def parse_seqtk_fqchk_output(seqtk_fqchk_output_path, quality_threshold):
    output = []
    with open(seqtk_fqchk_output_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            num_bases_and_avg_q = {}
            if row['position'] == 'ALL':
                percent_above_header = 'percent_bases_above_q' + str(quality_threshold)
                num_bases_and_avg_q['num_bases'] = int(row['num_bases'])
                num_bases_and_avg_q['average_q'] = float(row['average_q'])
                num_bases_and_avg_q[percent_above_header] = float(row[percent_above_header])
                output.append(num_bases_and_avg_q)
                
    return output      


def main(args):
    header = []
    with open(args.seqtk_fqchk_output,'r') as f:
        header = f.readline().strip().split(',')

    quality_threshold = 0
    for h in header:
        if h.startswith('percent_bases_above_q'):
            quality_threshold = int(h.split('percent_bases_above_q')[1])
            
    seqtk_fqchk_output = parse_seqtk_fqchk_output(args.seqtk_fqchk_output, quality_threshold)
    
    total_bases = sum([x['num_bases'] for x in seqtk_fqchk_output])

    overall_average_q = sum([x['average_q'] * x['num_bases'] for x in seqtk_fqchk_output]) / total_bases

    percent_above_header = 'percent_bases_above_q' + str(quality_threshold)

    overall_percent_above_threshold = sum([x[percent_above_header] * x['num_bases'] for x in seqtk_fqchk_output]) / total_bases

    print(','.join([
        'sample_id',
        'total_bases',
        'average_base_quality',
        percent_above_header,
    ]))
    
    print(','.join([
        args.sample_id,
        str(total_bases),
        str(round(overall_average_q, 3)),
        str(round(overall_percent_above_threshold, 3)),
    ]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('seqtk_fqchk_output')
    parser.add_argument('-s', '--sample-id')
    args = parser.parse_args()
    main(args)
