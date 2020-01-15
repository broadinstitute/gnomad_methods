#!/usr/bin/env bash -eu

DOCS_DIR=$(dirname "${BASH_SOURCE}")
cd $DOCS_DIR

rm -rf api_reference
./generate_api_reference.py

rm -rf html
python3 -m sphinx -W -b html . html
