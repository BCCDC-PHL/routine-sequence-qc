#!/bin/bash

set -eo pipefail

echo "Install ART .." >> artifacts/test_artifact.log

wget --quiet https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz && tar -xzf artbin*

mkdir -p /opt/art/bin

mv art_bin_*/art_illumina /opt/art/bin

echo "export PATH=/opt/art/bin:${PATH}" >> ~/.bashrc
