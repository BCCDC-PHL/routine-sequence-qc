#!/bin/bash

set -eo pipefail

export PATH=/opt/miniconda3/bin:${PATH}

echo "Download Reference Genomes..." >> artifacts/test_artifact.log

mkdir -p $PWD/.github/data/ref_genomes

pushd $PWD/.github/data/ref_genomes

while IFS=$'\t' read -r accession accession_version taxid gi ; do
    ncbi-acc-download --format fasta ${accession_version}
    sleep 5
done < <(tail -n+2 ../kraken2_db/taxonomy/nucl_gb.accession2taxid)

popd
