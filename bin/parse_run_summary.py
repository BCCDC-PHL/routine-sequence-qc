#!/usr/bin/python

import argparse
import collections
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
        '% Occupied': 'PercentOccupied',
        'Level': 'ReadNumber',
        'Intensity C1': 'IntensityCycle1',
    }

    level_to_read_number = {
        'Read 1': 1,
        'Read 2 (I)': 2,
        'Read 3 (I)': 3,
        'Read 4': 4,
    }

    with open(summary_path) as summary:
        for line in summary:
            if re.match("^Level", line):
                read_summary_headers = re.split("\s*,", line.rstrip())
                #read_summary_headers = [
                #    x.replace(" ", "") for x in read_summary_headers
                #]
                print(read_summary_headers)
                for idx, header in enumerate(read_summary_headers):
                    if header in replaced_fields:
                        read_summary_headers[idx] = replaced_fields[header]
                print(read_summary_headers)
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
        read_summary.append(read_summary_line_dict)
    
    return read_summary


def parse_lanes_by_read(summary_path):
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
                read_line = [record.strip() for record in line.strip().split(',')]
                read_line_dict = dict(zip(headers_input_order, read_line))
                read_line_dict['ReadNumber'] = read_number
            else:
                read_line_dict = {}

            if read_line_dict:
                for field in average_stdev_fields:
                    string_value = read_line_dict[field]
                    [average, stdev] = [float(value) for value in string_value.split(' +/- ')]
                    deviation_field = field + 'Deviation'
                    read_line_dict[field] = average
                    read_line_dict[deviation_field] = stdev

                for field, num_denom in slash_fields.items():
                    string_value = read_line_dict[field]
                    numerator_field = num_denom['numerator_field']
                    denominator_field = num_denom['denominator_field']
                    [numerator, denominator] = [float(value) for value in string_value.split(' / ')]
                    read_line_dict[numerator_field] = numerator
                    read_line_dict[denominator_field] = denominator
                    read_line_dict.pop(field, None)

                for field in float_fields:
                    read_line_dict[field] = float(read_line_dict[field])

                for field in int_fields:
                    read_line_dict[field] = int(read_line_dict[field])

                if read_line_dict['Surface'] == '-':
                    read_line_dict['Density'] = int(float(read_line_dict['Density'] * 1000))
                    read_line_dict['ClusterDensity'] = read_line_dict['Density']
                    read_line_dict['Reads'] = int(float(read_line_dict['Reads'] * 1000000))
                    read_line_dict['ReadsPf'] = int(float(read_line_dict['ReadsPf'] * 1000000))

                    read_line_dict.pop('OccupancyDeviation', None)
                    read_line_dict.pop('LegacyPhasingPrephasingRate', None)
                    read_line_dict.pop('Surface', None)

                    read_line_dict_ordered = collections.OrderedDict(sorted(read_line_dict.items(), key=lambda x: headers_output_order.index(x[0])))
                
                    lanes_by_read.append(read_line_dict_ordered)
    
    return lanes_by_read


def main(args):
    reads = parse_read_summary(args.summary)
    lanes_by_read = parse_lanes_by_read(args.summary)

    output = {
        'Reads': reads,
        'LanesByRead': lanes_by_read,
    }

    print(json.dumps(output, indent=2))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--summary')
    args = parser.parse_args()
    main(args)
