#!/usr/bin/env python

import argparse
import json
import re

from pprint import pprint


# https://stackoverflow.com/a/1176023
def camel_to_snake(name):
  name = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
  return re.sub('([a-z0-9])([A-Z])', r'\1_\2', name).lower()


def parse_header_section(path_to_sample_sheet):
  header_lines = []
  header = {}
  header['instrument_type'] = 'NextSeq2000'

  with open(path_to_sample_sheet, 'r') as f:
      for line in f:
          if line.strip().startswith('[Header]'):
              continue
          if line.strip().startswith('[Reads]') or line.strip().rstrip(',') == "":
              break
          else:
              header_lines.append(line.strip().rstrip(','))

  for line in header_lines:
      header_key = camel_to_snake(line.split(',')[0])

      if len(line.split(',')) > 1:
        header_value = line.split(',')[1]
      else:
        header_value = ""

      if header_key != "":
        header[header_key] = header_value
              
  return header


def parse_reads_section(path_to_sample_sheet):
  reads_lines = []
  reads = {}
  with open(path_to_sample_sheet, 'r') as f:
      for line in f:
          if line.strip().startswith('[Reads]'):
              break
      for line in f:
          if line.strip().startswith('[Sequencing_Settings]') or line.strip().rstrip(',') == "":
              break
          reads_lines.append(line.strip().rstrip(','))

  for line in reads_lines:
    reads_key = camel_to_snake(line.split(',')[0])
    if len(line.split(',')) > 1:
      reads_value = int(line.split(',')[1])
    else:
      reads_value = ""

    if reads_key != "":
      reads[reads_key] = reads_value

  return reads


def parse_sequencing_settings_section(path_to_sample_sheet):
  sequencing_settings_lines = []
  sequencing_settings = {}
  with open(path_to_sample_sheet, 'r') as f:
      for line in f:
          if line.strip().startswith('[Sequencing_Settings]'):
              break
      for line in f:
          if line.strip().startswith('[BCLConvert_Settings]') or line.strip().rstrip(',') == "":
              break
          sequencing_settings_lines.append(line.strip().rstrip(','))

  for line in sequencing_settings_lines:
    sequencing_settings_key = camel_to_snake(line.split(',')[0])

    if len(line.split(',')) > 1:
      sequencing_settings_value = line.split(',')[1]
    else:
      sequencing_settings_value = ""

    if sequencing_settings_key != "":
      sequencing_settings[sequencing_settings_key] = sequencing_settings_value
              
  return sequencing_settings


def parse_bclconvert_settings_section(path_to_sample_sheet):
  bclconvert_settings_lines = []
  bclconvert_settings = {}
  with open(path_to_sample_sheet, 'r') as f:
      for line in f:
          if line.strip().startswith('[BCLConvert_Settings]'):
              break
      for line in f:
          if line.strip().startswith('[BCLConvert_Data]') or line.strip().rstrip(',') == "":
              break
          bclconvert_settings_lines.append(line.strip().rstrip(','))

  for line in bclconvert_settings_lines:
    bclconvert_settings_key = camel_to_snake(line.split(',')[0])

    if len(line.split(',')) > 1:
      bclconvert_settings_value = line.split(',')[1]
    else:
      bclconvert_settings_value = ""

    if bclconvert_settings_key != "":
      bclconvert_settings[bclconvert_settings_key] = bclconvert_settings_value
              
  return bclconvert_settings


def parse_bclconvert_data_section(path_to_sample_sheet):
  bclconvert_data_lines = []
  bclconvert_data = []
  with open(path_to_sample_sheet, 'r') as f:
      for line in f:
          if line.strip().startswith('[BCLConvert_Data]'):
              break
      for line in f:
          if line.strip().startswith('[Cloud_Settings]') or line.strip().rstrip(',') == "":
              break
          bclconvert_data_lines.append(line.strip().rstrip(','))

  bclconvert_data_keys = [camel_to_snake(x) for x in bclconvert_data_lines[0].split(',')]
  for line in bclconvert_data_lines[1:]:
    d = {}
    values = line.split(',')
    for idx, key in enumerate(bclconvert_data_keys):
      d[key] = values[idx]
    bclconvert_data.append(d)

  return bclconvert_data

def parse_cloud_settings_section(path_to_sample_sheet):
  cloud_settings_lines = []
  cloud_settings = {}
  with open(path_to_sample_sheet, 'r') as f:
      for line in f:
          if line.strip().startswith('[Cloud_Settings]'):
              break
      for line in f:
          if line.strip().startswith('[BCLConvert_Settings]') or line.strip().rstrip(',') == "":
              break
          cloud_settings_lines.append(line.strip().rstrip(','))

  for line in cloud_settings_lines:
    cloud_settings_key = camel_to_snake(line.split(',')[0])

    if len(line.split(',')) > 1:
      cloud_settings_value = line.split(',')[1]
    else:
      cloud_settings_value = ""

    if cloud_settings_key != "":
      cloud_settings[cloud_settings_key] = cloud_settings_value
              
  return cloud_settings


def parse_cloud_data_section(path_to_sample_sheet):
  cloud_data_lines = []
  cloud_data = []
  with open(path_to_sample_sheet, 'r') as f:
      for line in f:
          if line.strip().startswith('[Cloud_Data]'):
              break
      for line in f:
          if line.strip().startswith('[Cloud_Settings]') or line.strip().rstrip(',') == "":
              break
          cloud_data_lines.append(line.strip().rstrip(','))

  if cloud_data_lines:
    cloud_data_keys = [camel_to_snake(x) for x in cloud_data_lines[0].split(',')]
    for line in cloud_data_lines[1:]:
      d = {}
      values = line.strip().split(',')
      if not all([x == '' for x in values]):
        for idx, key in enumerate(cloud_data_keys):
          try:
            d[key] = values[idx]
          except IndexError as e:
            d[key] = ""
          cloud_data.append(d)
  
  return cloud_data


def main(args):
  sample_sheet = {}
  header = parse_header_section(args.sample_sheet)
  reads = parse_reads_section(args.sample_sheet)
  sequencing_settings = parse_sequencing_settings_section(args.sample_sheet)
  bclconvert_settings = parse_bclconvert_settings_section(args.sample_sheet)
  bclconvert_data = parse_bclconvert_data_section(args.sample_sheet)
  cloud_settings = parse_cloud_settings_section(args.sample_sheet)
  cloud_data = parse_cloud_data_section(args.sample_sheet)
  sample_sheet['header'] = header
  sample_sheet['reads'] = reads
  sample_sheet['sequencing_settings'] = sequencing_settings
  sample_sheet['bclconvert_settings'] = bclconvert_settings
  sample_sheet['bclconvert_data'] = bclconvert_data
  if cloud_settings:
    sample_sheet['cloud_settings'] = cloud_settings
  if cloud_data:
    sample_sheet['cloud_data'] = cloud_data

  print(json.dumps(sample_sheet, indent=2))
  
  

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('sample_sheet')
  args = parser.parse_args()
  main(args)
