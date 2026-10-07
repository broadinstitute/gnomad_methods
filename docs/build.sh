#!/bin/sh

DOCS_DIR=$(dirname "$0")
cd "${DOCS_DIR}"

rm -rf api_reference
./generate_api_reference.py

# The knowledge tree lives at the repo root (readable on GitHub, not coupled
# to this build); stage a copy here so Sphinx can render it into the site.
rm -rf knowledge
cp -R ../knowledge knowledge

rm -rf html
python3 -m sphinx -W -b html . html
