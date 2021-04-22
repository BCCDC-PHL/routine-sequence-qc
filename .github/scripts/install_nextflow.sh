#!/bin/bash
set -eo pipefail

echo Install Nextflow .. >> artifacts/test_artifact.log

wget -qO- https://get.nextflow.io | bash > nextflow_install.out 2> nextflow_install.err

mkdir -p /opt/nextflow/bin

mv nextflow /opt/nextflow/bin

echo "export PATH=/opt/nextflow/bin:$PATH" >> ~/.bashrc

export PATH=/opt/nextflow/bin:$PATH

NXF_VER=20.10.0 nextflow run -C .github/data/config/nextflow.config hello

mv nextflow_install.out nextflow_install.err artifacts
