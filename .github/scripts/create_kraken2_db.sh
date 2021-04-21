#!/bin/bash

set -eo pipefail

export PATH=/opt/miniconda3/bin:$PATH

pushd ${PWD}/.github/data && mkdir -p kraken2_db

kraken2-build --download-taxonomy --db kraken2_db

for file in ref_genomes/*.fa; do
    kraken2-build --add-to-library ${file} --db kraken2_db
done

kraken2-build --build --db kraken2_db > kraken_build.log

popd

cp ${PWD}/.github/data/kraken_build.log artifacts
