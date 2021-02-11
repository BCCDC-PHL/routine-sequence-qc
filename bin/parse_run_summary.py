#!/usr/bin/python

import argparse
import sys
import re
import json

def parse_read_summary(summary_path):
    read_summary_headers = []
    read_summary_lines = []
    # Basic approach to parsing text between two specific lines
    # described here: https://stackoverflow.com/a/7559542/780188
    with open(summary_path) as summary:
        for line in summary:
            if re.match("^Level", line):
                read_summary_headers = re.split("\s*,", line.rstrip())
                read_summary_headers = [
                    x.lower().replace(" ", "_") for x in read_summary_headers
                ]
                read_summary_headers = [
                    x.replace("%>=q30", "percent_greater_than_q30") for x in read_summary_headers
                ]
                break
        for line in summary:
            if re.match("^Total", line):
                read_summary_lines.append(re.split("\s*,", line.rstrip()))
                break
            read_summary_lines.append(re.split("\s*,", line.rstrip()))

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
        'legacy_pasing_prephasing_rate',
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
    lines_by_read = {
        'read_1': [],
        'read_i1': [],
        'read_i2': [],
        'read_2': [],
    }
    with open(summary_path) as summary:
        current_read = None
        for line in summary:
            if re.match("^Read 1$", line):
                current_read = 'read_1'
            elif re.match("^Read 2 \(I\)$", line):
                current_read = 'read_i1'
            elif re.match("^Read 3 \(I\)$", line):
                current_read = 'read_i2'
            elif re.match("^Read 4$", line):
                current_read = 'read_2'
            elif re.match("^Read 4$", line):
                current_read = 'read_2'
            elif re.match("^Extracted", line) or re.match("^Called", line) or re.match("^Scored", line):
                current_read = None
            elif current_read and not re.match("^Lane", line):
                read_line = [record.strip() for record in line.strip().split(',')]
                read_line_dict = dict(zip(headers, read_line))
                lines_by_read[current_read].append(read_line_dict)
            else:
                pass
    
    return lines_by_read

def main(args):
    read_summary = parse_read_summary(args.summary)
    read_summary_detail = parse_read_summary_detail(args.summary)

    # print(json.dumps(read_summary))
    print(json.dumps(read_summary_detail))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--summary')
    args = parser.parse_args()
    main(args)
