#!/bin/bash


ROLE=$(/usr/share/google/get_metadata_value attributes/dataproc-role)
if [[ "${ROLE}" == 'Master' ]]; then

    easy_install pip
    pip install decorator

    mkdir -p /home/hail/

    git clone git@github.com:macarthur-lab/gnomad_hail.git /home/hail/
    cd /home/hail/gnomad_hail
    chmod +x ./init_scripts/gnomad-init.sh ./init_scripts/sparklyr-init.sh
    ./init_scripts/gnomad-init.sh > ./init_scripts/gnomad_startup.log 2>&1 &
    ./init_scripts/sparklyr-init.sh > ./init_scripts/sparklyr_startup.log 2>&1 &

    python -m unittest discover &> tests.log

fi