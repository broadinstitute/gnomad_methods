#!/bin/bash


ROLE=$(/usr/share/google/get_metadata_value attributes/dataproc-role)
if [[ "${ROLE}" == 'Master' ]]; then

    easy_install pip
    pip install decorator
    gsutil cp gs://gnomad-public/tools/inits/gnomad-init.sh .
    chmod +x gnomad-init.sh
    ./gnomad-init.sh > gnomad_startup.log 2>&1 &

    gsutil cp gs://gnomad-public/tools/inits/sparklyr-init.sh .
    chmod +x sparklyr-init.sh
    ./sparklyr-init.sh > sparklyr_startup.log 2>&1 &
fi