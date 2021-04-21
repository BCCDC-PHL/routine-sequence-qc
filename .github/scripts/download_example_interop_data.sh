#!/bin/bash

set -eo pipefail

pushd ${PWD}/.github/data

# Only publicly-available InterOp data I've found...
wget http://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz

tar -xzf cellranger-tiny-bcl-1.2.0.tar.gz && rm cellranger-tiny-bcl-1.2.0.tar.gz

mv cellranger-tiny-bcl-1.2.0/InterOp mock_runs/210101_M00000_0000_000000000-A1B2C
mv cellranger-tiny-bcl-1.2.0/RunInfo.xml mock_runs/210101_M00000_0000_000000000-A1B2C
mv cellranger-tiny-bcl-1.2.0/runParameters.xml mock_runs/210101_M00000_0000_000000000-A1B2C

rm -r cellranger-tiny-bcl-1.2.0

popd
