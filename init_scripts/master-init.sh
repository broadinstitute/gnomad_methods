#!/bin/bash


ROLE=$(/usr/share/google/get_metadata_value attributes/dataproc-role)
if [[ "${ROLE}" == 'Master' ]]; then

    easy_install pip
    pip install decorator

    mkdir -p /home/hail/

    cd /home/hail
    git clone https://github.com/macarthur-lab/gnomad_hail.git
    cd gnomad_hail/init_scripts
    chmod +x gnomad-init.sh sparklyr-init.sh
    ./gnomad-init.sh > gnomad_startup.log 2>&1 &
    ./sparklyr-init.sh > sparklyr_startup.log 2>&1 &

fi