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

# RunInfo.xml doesn't match InterOp data for some reason, so
# edit RunInfo.xml to match
cat mock_runs/210101_M00000_0000_000000000-A1B2C/RunInfo.xml | \
    sed 's/LaneCount="1"/LaneCount="2"/' | \
    sed 's/SurfaceCount="1"/SurfaceCount="2"/' | \
    sed 's/SwathCount="1"/SwathCount="2"/' | \
    sed 's/TileCount="1"/TileCount="32"/' \
	> RunInfo.edited.xml

mv RunInfo.edited.xml mock_runs/210101_M00000_0000_000000000-A1B2C/RunInfo.xml

popd
