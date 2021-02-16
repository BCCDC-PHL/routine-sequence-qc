#!/usr/bin/python

import argparse
import sys
import re
import json

def parse_read_summary(summary_path):
    read_summary_headers = []
    read_summary_lines = []

    replaced_fields = {'%>=q30': 'percent_greater_than_q30',
                       '%_occupied': 'percent_occupied'}

    with open(summary_path) as summary:
        for line in summary:
            if re.match("^Level", line):
                read_summary_headers = re.split("\s*,", line.rstrip())
                read_summary_headers = [
                    x.lower().replace(" ", "_") for x in read_summary_headers
                ]
                for idx, header in enumerate(read_summary_headers):
                    if header in replaced_fields:
                        read_summary_headers[idx] = replaced_fields[header]
                
                break
        for line in summary:
            if re.match("^Total", line):
                read_summary_lines.append(re.split(",", line.rstrip()))
                break
            else:
                read_summary_lines.append(re.split(",", line.rstrip()))

    read_summary = []
    for line in read_summary_lines:
        read_summary_line_dict = {}
        for idx, header in enumerate(read_summary_headers):
            if header == 'level':
                read_summary_line_dict[header] = line[idx]
            elif header == 'intensity_c1':
                read_summary_line_dict[header] = int(line[idx])
            else:
                read_summary_line_dict[header] = float(line[idx])
        read_summary.append(read_summary_line_dict)
    
    return read_summary


def parse_read_summary_detail(summary_path):
    headers = [
        'lane',
        'surface',
        'tiles',
        'density',
        'clusters_passing_filter',
        'legacy_phasing_prephasing_rate',
        'phasing_slope_offset',
        'prephasing_slope_offset',
        'reads',
        'reads_passing_filter',
        'percent_greater_than_q30',
        'yield',
        'cycles_error',
        'aligned',
        'error',
        'error_35',
        'error_75',
        'error_100',
        'percent_occupied',
        'intensity_at_cycle_1',
    ]
    average_stdev_fields = [
        'aligned',
        'clusters_passing_filter',
        'density',
        'error',
        'error_100',
        'error_75',
        'error_35',
        'intensity_at_cycle_1',
        'percent_occupied',
    ]
    slash_fields = { 'legacy_phasing_prephasing_rate': {'numerator_field': 'legacy_phasing_rate',
                                                        'denominator_field': 'legacy_prephasing_rate'},
                     'phasing_slope_offset': {'numerator_field': 'phasing_slope',
                                              'denominator_field': 'phasing_offset'},
                     'prephasing_slope_offset': {'numerator_field': 'prephasing_slope',
                                                 'denominator_field': 'prephasing_offset'},
    }
    float_fields = [
        'percent_greater_than_q30',
        'reads',
        'reads_passing_filter',
        'yield',
    ]
    
    lines_by_read = {
        'read_1': [],
        'read_i1': [],
        'read_i2': [],
        'read_2': [],
    }
    with open(summary_path) as summary:
        current_read = None
        for line in summary:
            if re.match("^Read 1\n$", line):
                current_read = 'read_1'
            elif re.match("^Read 2 \(I\)\n$", line):
                current_read = 'read_i1'
            elif re.match("^Read 3 \(I\)\n$", line):
                current_read = 'read_i2'
            elif re.match("^Read 4$\n", line):
                current_read = 'read_2'
            elif re.match("^Extracted", line) or re.match("^Called", line) or re.match("^Scored", line):
                current_read = None
            elif current_read and not re.match("^Lane", line):
                read_line = [record.strip() for record in line.strip().split(',')]
                read_line_dict = dict(zip(headers, read_line))
                lines_by_read[current_read].append(read_line_dict)
            else:
                pass

        for field in average_stdev_fields:
            string_value = read_line_dict[field]
            [average, stdev] = [float(value) for value in string_value.split(' +/- ')]
            read_line_dict[field] = { 'average': average,
                                      'stdev': stdev }

        for field, num_denom in slash_fields.items():
            string_value = read_line_dict[field]
            numerator_field = num_denom['numerator_field']
            denominator_field = num_denom['denominator_field']
            [numerator, denominator] = [float(value) for value in string_value.split(' / ')]
            read_line_dict[numerator_field] = numerator
            read_line_dict[denominator_field] = denominator
            read_line_dict.pop(field, None)
    
    return lines_by_read


def main(args):
    read_summary = parse_read_summary(args.summary)
    read_summary_detail = parse_read_summary_detail(args.summary)

    output = {'read_summary': read_summary,
              'read_details': read_summary_detail}
    # print(json.dumps(read_summary))
    print(json.dumps(output))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--summary')
    args = parser.parse_args()
    main(args)
