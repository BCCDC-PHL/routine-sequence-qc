#!/bin/bash

set -eo pipefail

pushd ${PWD}/.github/data/mock_runs/210101_M00000_0000_000000000-A1B2C/Data/Intensities/BaseCalls/

art_illumina --seqSys MSv3 --paired -i ../../../../../ref_genomes/NC_002695.2.fa --fcov 5 --mflen 500 --sdev 25 --len 250 --noALN -o test-01_R

art_illumina --seqSys MSv3 --paired -i ../../../../../ref_genomes/NC_016845.1.fa --fcov 5 --mflen 500 --sdev 25 --len 250 --noALN -o test-02_R

art_illumina --seqSys MSv3 --paired -i ../../../../../ref_genomes/NZ_CP033744.1.fa --fcov 5 --mflen 500 --sdev 25 --len 250 --noALN -o test-03_R

art_illumina --seqSys MSv3 --paired -i ../../../../../ref_genomes/NC_003197.2.fa --fcov 5 --mflen 500 --sdev 25 --len 250 --noALN -o negative-control_R

rename fq fastq *.fq

gzip *.fastq

popd
