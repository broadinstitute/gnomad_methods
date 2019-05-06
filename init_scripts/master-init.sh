#!/bin/bash

set -x


# Installing numpy on workers for pyspark purposes
/opt/conda/bin/pip install numpy

ROLE=$(/usr/share/google/get_metadata_value attributes/dataproc-role)
if [[ "${ROLE}" == 'Master' ]]; then

    mkdir -p /home/hail/

    cd /home/hail
    git clone https://github.com/macarthur-lab/gnomad_hail.git
    cd gnomad_hail/init_scripts
    chmod +x gnomad-init.sh sparklyr-init.sh

    # This is here so as not have 2 apt-get processes fighting for a lock
    apt-get install -y tmux liblz4-dev

    ./gnomad-init.sh > gnomad_startup.log 2>&1 &
    ./sparklyr-init.sh > sparklyr_startup.log 2>&1 &

fi
