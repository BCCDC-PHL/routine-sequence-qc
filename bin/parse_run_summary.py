#!/usr/bin/python

import argparse
import collections
import math
import sys
import re
import json

def parse_read_summary(summary_path):
    read_summary_headers = []
    read_summary_lines = []

    replaced_fields = {
        '%>=Q30': 'PercentGtQ30',
        'Projected Yield': 'ProjectedTotalYield',
        'Yield': 'YieldTotal',
        'Error Rate': 'ErrorRate',
        'Aligned': 'PercentAligned',
        '% Occupied': 'Occupancy',
        'Level': 'ReadNumber',
        'Intensity C1': 'IntensityCycle1',
    }

    level_to_read_number = {
        'Read 1': 1,
        'Read 2 (I)': 2,
        'Read 3 (I)': 3,
        'Read 4': 4,
    }

    headers_output_order = [
        'ReadNumber',
        'IsIndexed',
        'TotalCycles',
        'YieldTotal',
        'ProjectedTotalYield',
        'PercentAligned',
        'ErrorRate',
        'IntensityCycle1',
        'PercentGtQ30',
    ]

    with open(summary_path) as summary:
        for line in summary:
            if re.match("^Level", line):
                read_summary_headers = re.split("\s*,", line.rstrip())
                
                for idx, header in enumerate(read_summary_headers):
                    if header in replaced_fields:
                        read_summary_headers[idx] = replaced_fields[header]
                break
        for line in summary:
            if re.match("^Total", line) or re.match("^Non-indexed", line):
                # read_summary_lines.append(re.split(",", line.rstrip()))
                break
            else:
                read_summary_lines.append(re.split(",", line.rstrip()))

    read_summary = []
    for line in read_summary_lines:
        read_summary_line_dict = {}
        for idx, header in enumerate(read_summary_headers):
            if header == 'ReadNumber':
                read_summary_line_dict[header] = level_to_read_number[line[idx]]
                if re.search("(I)", line[idx]):
                    read_summary_line_dict['IsIndexed'] = True
                else:
                    read_summary_line_dict['IsIndexed'] = False
            elif header == 'IntensityC1':
                read_summary_line_dict[header] = int(line[idx])
            elif header == 'ErrorRate':
                if line[idx] == 'nan':
                    read_summary_line_dict[header] = 0
                else:
                    read_summary_line_dict[header] = float(line[idx])
            else:
                read_summary_line_dict[header] = float(line[idx])

        read_summary_line_dict.pop('Occupancy', None)
        read_summary_line_dict_ordered = collections.OrderedDict(sorted(read_summary_line_dict.items(), key=lambda x: headers_output_order.index(x[0])))

        read_summary.append(read_summary_line_dict_ordered)
    
    return read_summary


def parse_read_line(read_line, read_number):
    parsed_read_line = {}

    headers_input_order = [
        'LaneNumber',
        'Surface',
        'TileCount',
        'Density',
        'PercentPf',
        'LegacyPhasingPrephasingRate',
        'Phasing',
        'Prephasing',
        'Reads',
        'ReadsPf',
        'PercentGtQ30',
        'Yield',
        'CyclesError',
        'PercentAligned',
        'ErrorRate',
        'ErrorRate35',
        'ErrorRate75',
        'ErrorRate100',
        'Occupancy',
        'IntensityCycle1',
    ]

    headers_output_order = [
        'ReadNumber',
        'LaneNumber',
        'Surface',
        'TileCount',
        'Density',
        'DensityDeviation',
        'PercentPf',
        'PercentPfDeviation',
        'Reads',
        'ReadsPf',
        'PercentGtQ30',
        'Yield',
        'CyclesError',
        'PercentAligned',
        'PercentAlignedDeviation',
        'ErrorRate',
        'ErrorRateDeviation',
        'ErrorRate35',
        'ErrorRate35Deviation',
        'ErrorRate75',
        'ErrorRate75Deviation',
        'ErrorRate100',
        'ErrorRate100Deviation',
        'IntensityCycle1',
        'IntensityCycle1Deviation',
        'PhasingSlope',
        'PhasingOffset',
        'PrePhasingSlope',
        'PrePhasingOffset',
        'ClusterDensity',
        'Occupancy',
    ]

    average_stdev_fields = [
        'PercentAligned',
        'PercentPf',
        'Density',
        'ErrorRate',
        'ErrorRate35',
        'ErrorRate75',
        'ErrorRate100',
        'IntensityCycle1',
        'Occupancy',
    ]

    slash_fields = { 
        'Phasing': {
            'numerator_field': 'PhasingSlope',
            'denominator_field': 'PhasingOffset'
        },
        'Prephasing': {
            'numerator_field': 'PrePhasingSlope',
            'denominator_field': 'PrePhasingOffset'
        },
    }

    float_fields = [
        'PercentGtQ30',
        'Reads',
        'ReadsPf',
        'Yield',
    ]

    int_fields = [
        'ReadNumber',
        'LaneNumber',
        'TileCount',
    ]

    read_line_list = [record.strip() for record in read_line.strip().split(',')]
    parsed_read_line = dict(zip(headers_input_order, read_line_list))
    parsed_read_line['ReadNumber'] = read_number

    for field in average_stdev_fields:
        string_value = parsed_read_line[field]
        [average, stdev] = [float(value) for value in string_value.split(' +/- ')]
        deviation_field = field + 'Deviation'
        parsed_read_line[field] = average
        parsed_read_line[deviation_field] = stdev

    for field, num_denom in slash_fields.items():
        string_value = parsed_read_line[field]
        numerator_field = num_denom['numerator_field']
        denominator_field = num_denom['denominator_field']
        [numerator, denominator] = [float(value) for value in string_value.split(' / ')]
        parsed_read_line[numerator_field] = numerator
        parsed_read_line[denominator_field] = denominator
        parsed_read_line.pop(field, None)

    for field in float_fields:
        parsed_read_line[field] = float(parsed_read_line[field])

    for field in int_fields:
        parsed_read_line[field] = int(parsed_read_line[field])

    parsed_read_line['Density'] = int(float(parsed_read_line['Density'] * 1000))
    parsed_read_line['ClusterDensity'] = parsed_read_line['Density']
    parsed_read_line['Reads'] = int(float(parsed_read_line['Reads'] * 1000000))
    parsed_read_line['ReadsPf'] = int(float(parsed_read_line['ReadsPf'] * 1000000))

    parsed_read_line.pop('OccupancyDeviation', None)
    parsed_read_line.pop('LegacyPhasingPrephasingRate', None)

    for k, v in parsed_read_line.items():
        if type(v) is float and math.isnan(v):
            parsed_read_line[k] = 0

    parsed_read_line_ordered = collections.OrderedDict(sorted(parsed_read_line.items(), key=lambda x: headers_output_order.index(x[0])))
    
    return parsed_read_line_ordered


def parse_lanes_by_read(summary_path):    
    lanes_by_read = []

    with open(summary_path) as summary:
        read_number = None
        for line in summary:
            if re.match("^Read 1\n$", line):
                read_number = 1
            elif re.match("^Read 2 \(I\)\n$", line):
                read_number = 2
            elif re.match("^Read 3 \(I\)\n$", line):
                read_number = 3
            elif re.match("^Read 4$\n", line):
                read_number = 4
            elif re.match("^Extracted", line) or re.match("^Called", line) or re.match("^Scored", line):
                read_number = None
            else:
                pass

            if read_number and not re.match("^Lane", line) and not re.match("^Read", line):
                parsed_read_line = parse_read_line(line, read_number)
                parsed_read_line['ReadNumber'] = read_number
            else:
                parsed_read_line = {}

            if parsed_read_line and parsed_read_line['Surface'] == '-':
                parsed_read_line.pop('Surface', None)
                lanes_by_read.append(parsed_read_line)
    
    return lanes_by_read




def main(args):
    reads = parse_read_summary(args.summary)
    lanes_by_read = parse_lanes_by_read(args.summary)

    output = collections.OrderedDict()
    output['Reads'] = reads
    output['LanesByRead'] = lanes_by_read

    print(json.dumps(output, indent=2))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--summary')
    args = parser.parse_args()
    main(args)
