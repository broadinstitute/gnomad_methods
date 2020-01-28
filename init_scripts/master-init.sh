#!/bin/bash

set -x


# Installing numpy on workers for pyspark purposes
pip install numpy

ROLE=$(/usr/share/google/get_metadata_value attributes/dataproc-role)
if [[ "${ROLE}" == 'Master' ]]; then

    git clone https://github.com/macarthur-lab/gnomad_hail.git /home/gnomad_hail

    mkdir -p /home/hail/
    ln -s ln -s /home/gnomad_hail/gnomad_hail /home/hail

    cd /home/gnomad_hail/init_scripts
    chmod +x gnomad-init.sh sparklyr-init.sh

    # This is here so as not have 2 apt-get processes fighting for a lock
    apt-get install -y tmux liblz4-dev

    ./gnomad-init.sh > gnomad_startup.log 2>&1 &
    ./sparklyr-init.sh > sparklyr_startup.log 2>&1 &

fi
