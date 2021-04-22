#!/bin/bash
set -eo pipefail

echo "Install System Dependencies .." >> artifacts/test_artifact.log

sudo apt-get update --fix-missing

sudo apt-get install -y \
     wget \
     bzip2 \
     ca-certificates \
     libglib2.0-0 \
     libxext6 \
     libsm6 libxrender1 \
     git \
     mercurial \
     subversion \
     rename
