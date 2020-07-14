#!/bin/sh

DOCS_DIR=$(dirname "$0")
cd "${DOCS_DIR}"

rm -rf api_reference
./generate_api_reference.py

./copy_changelog.py

rm -rf html
python3 -m sphinx -W -b html . html
