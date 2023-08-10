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
            try:
                row['fraction_total_reads'] = float(row['fraction_total_reads'])
            except ValueError as e:
                row['fraction_total_reads'] = None
            bracken_report_lines.append(row)

    return bracken_report_lines
        

def main(args):
    bracken_report = parse_bracken_report(args.bracken_report)

    fraction_unclassified_reads = 0.0
    unclassified_reads_records = list(filter(lambda x: x['taxonomy_id'] == '0', bracken_report))
    if len(unclassified_reads_records) > 0:
        unclassified_reads_record = unclassified_reads_records[0]
        fraction_unclassified_reads = unclassified_reads_record['fraction_total_reads']

    bracken_report_without_unclassified = list(filter(lambda x: x['taxonomy_id'] != '0', bracken_report))
    bracken_report_sorted = sorted(bracken_report_without_unclassified, key=lambda k: k['fraction_total_reads'], reverse=True) 
    
    output_fields = ['sample_id', 'taxonomy_level']
    output_line = {
        'sample_id': args.sample_id,
        'taxonomy_level': args.taxonomy_level,
    }

    for n in range(args.top_n):
        num = str(n + 1)
        name_field = 'abundance_' + num + '_name'
        try:
            output_line[name_field] = bracken_report_sorted[n]['name']
        except IndexError as e:
            output_line[name_field] = None
        output_fields.append(name_field)

        taxid_field = 'abundance_' + num + '_taxonomy_id'
        try:
            output_line[taxid_field] = bracken_report_sorted[n]['taxonomy_id']
        except IndexError as e:
            output_line[taxid_field] = None
        output_fields.append(taxid_field)

        fraction_total_reads_field = 'abundance_' + num + '_fraction_total_reads'
        try:
            output_line[fraction_total_reads_field] = round(bracken_report_sorted[n]['fraction_total_reads'], 6)
        except IndexError as e:
            output_line[fraction_total_reads_field] = 0.0
        output_fields.append(fraction_total_reads_field)

    unclassified_fields = [
        'abundance_unclassified_name',
        'abundance_unclassified_taxonomy_id',
        'abundance_unclassified_fraction_total_reads',
    ]
    output_fields += unclassified_fields
    output_line['abundance_unclassified_name'] = 'unclassified'
    output_line['abundance_unclassified_taxonomy_id'] = '0'
    output_line['abundance_unclassified_fraction_total_reads'] = fraction_unclassified_reads

    writer = csv.DictWriter(sys.stdout, fieldnames=output_fields, dialect='unix', quoting=csv.QUOTE_MINIMAL)
    writer.writeheader()
    writer.writerow(output_line)
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('bracken_report')
    parser.add_argument('-s', '--sample-id')
    parser.add_argument('-l', '--taxonomy-level')
    parser.add_argument('-n', '--top-n', type=int)
    args = parser.parse_args()
    main(args)
