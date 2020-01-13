#!/usr/bin/env bash -eu

DOCS_DIR=$(dirname "${BASH_SOURCE}")
cd $DOCS_DIR

python3 -m sphinx -W -b html . html
