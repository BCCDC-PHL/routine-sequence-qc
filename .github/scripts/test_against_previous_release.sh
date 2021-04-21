#!/bin/bash

set -eo pipefail

export PATH=/opt/miniconda3/bin:$PATH
export PATH=/opt/nextflow/bin:$PATH

# write test log as github Action artifact
echo "Nextflow run current PR..." >> artifacts/test_artifact.log
NXF_VER=20.10.0 nextflow -quiet run ./main.nf \
       -profile conda \
       --cache ~/.conda/envs \
       --kraken2_db $PWD/.github/data/kraken2_db \
       --bracken_db $PWD/.github/data/kraken2_db \
       --run_dir $PWD/.github/data/mock_runs/210101_M00000_0000_000000000-A1B2C \
       --outdir pr_output

cp .nextflow.log artifacts/

# run tests against previous previous_release to compare outputs 
git clone https://github.com/BCCDC-PHL/routine-sequence-qc.git previous_release 
cd previous_release
git checkout 26220ef1217229beb73393e74c56a57ea90150bf

echo "Nextflow run previous release..." >> ../artifacts/test_artifact.log
NXF_VER=20.10.0 nextflow -quiet run ./main.nf \
       -profile conda \
       --cache ~/.conda/envs \
       --directory $PWD/../.github/data/fastqs/ \
       --ref $PWD/../.github/data/refs/MN908947.3/MN908947.3.fa \
       --bed $PWD/../.github/data/primer_schemes/nCoV-2019_Freed_1200bp.bed \
       --primer_pairs_tsv $PWD/../.github/data/primer_schemes/nCoV-2019_Freed_1200bp_primer_pairs.tsv \
       --gff $PWD/../.github/data/refs/MN908947.3.gff \
       --composite_ref $PWD/../.github/data/refs/mock_composite_ref/mock_composite_ref.fa \
       --illumina \
       --prefix test

cp .nextflow.log ../artifacts/previous_release.nextflow.log

cd ..

# exclude files from comparison
# and list differences
echo "Compare ouputs of current PR vs those of previous release.." >> artifacts/test_artifact.log
find results ./previous_release/results \
     -name "*.fq.gz" \
     -o -name "*.bam" \
     -o -name "*.bam.bai" \
     -o -name "*.vcf" \
    | xargs rm -rf
if ! git diff --stat --no-index results ./previous_release/results > diffs.txt ; then
  echo "test failed: differences found between PR and previous release" >> artifacts/test_artifact.log
  echo "see diffs.txt" >> artifacts/test_artifact.log 
  cp diffs.txt artifacts/  
  exit 1
else
  echo "no differences found between PR and previous release" >> artifacts/test_artifact.log
fi

# clean-up for following tests
rm -rf previous_release && rm -rf results && rm -rf work && rm -rf .nextflow*
