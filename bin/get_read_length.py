#!/usr/bin/env python

import argparse
import json


def main(args):
  with open(args.sample_sheet_json, 'r') as f:
    sample_sheet = json.load(f)

  if sample_sheet['header']['instrument_type'] == "MiSeq":
    actual_read_length = sample_sheet['reads'][0]
  elif sample_sheet['header']['instrument_type'] == "NextSeq2000":
    actual_read_length = sample_sheet['reads']['read1_cycles']

  if actual_read_length < 60:
      read_length = 50
  elif actual_read_length < 90:
      read_length = 75
  elif actual_read_length < 125:
    read_length = 100
  elif actual_read_length < 175:
    read_length = 150
  elif actual_read_length < 225:
    read_length = 200
  elif actual_read_length < 275:
    read_length = 250
  else:
    read_length = 300
  print(read_length)


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('sample_sheet_json')
  args = parser.parse_args()
  main(args)
