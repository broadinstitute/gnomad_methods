#!/usr/bin/env bash

cat <<EOT >> /config.yaml
default:
  spark.kryoserializer.buffer.max: 1g
  spark.yarn.executor.memoryOverhead: 1861
  spark.executor.memory: 18619m
  spark.yarn.am.memoryOverhead: 1861
  spark.driver.memory: 45g
  spark.executor.cores: 4
  spark.yarn.am.memory: 18619m
  spark.driver.maxResultSize: 30g
  spark.task.maxFailures: 20
EOT

apt-get install -y gdebi-core libcurl4-openssl-dev libssl-dev libxml2-dev
export RSTUDIO_DEB=rstudio-server-1.0.143-amd64.deb
wget https://download2.rstudio.org/${RSTUDIO_DEB}
gdebi -n ${RSTUDIO_DEB}
mkdir /opt/spark
ln -s /usr/lib/spark /opt/spark/spark-2.0.2-bin-hadoop2.7

# Common R packages, tiered for faster startup
R --vanilla -e "install.packages(c('sparklyr', 'dplyr'), repos='https://cran.rstudio.com')"
R --vanilla -e "install.packages(c('magrittr', 'ggplot2', 'slackr', 'ggrepel'), repos='https://cran.rstudio.com')"
R --vanilla -e "install.packages(c('plyr', 'shiny', 'plotly'), repos='https://cran.rstudio.com')"
R --vanilla -e "install.packages(c('DT', 'tidyverse', 'broom', 'randomForest', 'ROCR', 'shinythemes', 'devtools'), repos='https://cran.rstudio.com')"
