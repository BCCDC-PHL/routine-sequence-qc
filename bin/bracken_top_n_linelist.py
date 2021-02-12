#!/usr/bin/env python

import argparse
import csv
import sys
import re
import json


def parse_bracken_report(bracken_report_path):
    bracken_report_lines = []
    with open(bracken_report_path, 'r') as f:
        reader = csv.DictReader(f, dialect='excel-tab')
        for row in reader:
            bracken_report_lines.append(row)

    return bracken_report_lines
        

def main(args):
    bracken_report = parse_bracken_report(args.bracken_report)
    
    # print(json.dumps(bracken_report))
    output_fields = ['sample_id', 'taxonomy_level']
    output_line = {
        'sample_id': args.sample_id,
        'taxonomy_level': bracken_report[0]['taxonomy_lvl']
    }
    
    for n in range(args.top_n):
        num = str(n + 1)
        name_field = 'abundance_' + num + '_name'
        output_line[name_field] = bracken_report[n]['name']
        output_fields.append(name_field)
        taxonomy_id_field = 'abundance_' + num + '_taxonomy_id'
        output_line[taxonomy_id_field] = bracken_report[n]['taxonomy_id']
        output_fields.append(taxonomy_id_field)
        estimated_reads_field = 'abundance_' + num + '_estimated_reads'
        output_line[estimated_reads_field] = bracken_report[n]['new_est_reads']
        output_fields.append(estimated_reads_field)
        fraction_total_reads_field = 'abundance_' + num + '_fraction_total_reads'
        output_line[fraction_total_reads_field] = bracken_report[n]['fraction_total_reads']
        output_fields.append(fraction_total_reads_field)
        

    writer = csv.DictWriter(sys.stdout, fieldnames=output_fields)
    writer.writeheader()
    writer.writerow(output_line)
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('bracken_report')
    parser.add_argument('-s', '--sample-id')
    parser.add_argument('-n', '--top-n', type=int)
    args = parser.parse_args()
    main(args)
