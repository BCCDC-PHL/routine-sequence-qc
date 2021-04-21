#!/bin/bash

set -eo pipefail

export PATH=/opt/miniconda3/bin:$PATH

pushd ${PWD}/.github/data

bracken-build -d kraken2_db -l 250 > bracken_build.log

kraken2-build --clean --db kraken2_db

popd

cp ${PWD}/.github/data/bracken_build.log artifacts
