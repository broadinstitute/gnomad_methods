#!/bin/bash


ROLE=$(/usr/share/google/get_metadata_value attributes/dataproc-role)
if [[ "${ROLE}" == 'Master' ]]; then

    mkdir -p /home/hail/

    cd /home/hail
    git clone https://github.com/macarthur-lab/gnomad_hail.git
    cd gnomad_hail
    git checkout hail-0.2-python
    cd init_scripts
    chmod +x gnomad-init-devel-spark2.2.0.sh sparklyr-init.sh

    # This is here so as not have 2 apt-get processes fighting for a lock
    apt-get install -y ipython tmux

    ./gnomad-init-devel.sh > gnomad_startup.log 2>&1 &
    ./sparklyr-init.sh > sparklyr_startup.log 2>&1 &

fi