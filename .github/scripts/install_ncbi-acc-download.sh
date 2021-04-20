#!/bin/bash

set -eo pipefail

echo "Install NCBI Genome Download tool using conda.." >> artifacts/test_artifact.log

conda install -y ncbi-acc-download
