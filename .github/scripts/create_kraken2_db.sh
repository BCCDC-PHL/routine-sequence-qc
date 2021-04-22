#!/bin/bash

set -eo pipefail

export PATH=/opt/miniconda3/bin:$PATH

pushd ${PWD}/.github/data/kraken2_db/taxonomy

rsync --no-motd rsync://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz .

tar -xzf taxdump.tar.gz && rm taxdump.tar.gz

pushd ../..

for file in ref_genomes/*.fa; do
    kraken2-build --add-to-library ${file} --db kraken2_db
done

kraken2-build --build --db kraken2_db > kraken_build.log

popd && popd

cp ${PWD}/.github/data/kraken_build.log artifacts
