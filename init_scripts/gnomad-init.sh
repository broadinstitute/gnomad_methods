#!/usr/bin/env bash

set -x

pip install --upgrade Cython
PACKAGES="slackclient==2.0.0 sklearn tabulate scipy statsmodels ggplot hdbscan websocket-client scikit-learn joblib"
pip install --upgrade --no-deps $PACKAGES

#export HAIL_VERSION=devel
#export SPARK_VERSION=2.2.0
#export SPARK_HOME=/usr/lib/spark
#export HAIL_HOME=/hadoop_gcs_connector_metadata_cache/hail
#export HAIL_HASH=$(gsutil cat gs://hail-common/builds/${HAIL_VERSION}/latest-hash-spark-${SPARK_VERSION}.txt)
#export HAIL_JAR=hail-${HAIL_VERSION}-${HAIL_HASH}-Spark-${SPARK_VERSION}.jar
#export HAIL_PYTHON_ZIP=hail-${HAIL_VERSION}-${HAIL_HASH}.zip
#mkdir $HAIL_HOME
#gsutil cp gs://hail-common/builds/${HAIL_VERSION}/jars/${HAIL_JAR} gs://hail-common/builds/${HAIL_VERSION}/python/${HAIL_PYTHON_ZIP} $HAIL_HOME
#
#mkdir -p /home/hail/conf-ipython
#cp /etc/spark/conf/spark-defaults.conf /etc/spark/conf/spark-env.sh /home/hail/conf-ipython/
#cat <<EOT >> /home/hail/conf-ipython/spark-defaults.conf
#spark.files=${HAIL_HOME}/${HAIL_JAR}
#spark.submit.pyFiles=${HAIL_HOME}/${HAIL_PYTHON_ZIP}
#spark.driver.extraClassPath=${HAIL_HOME}/${HAIL_JAR}
#spark.executor.extraClassPath=${HAIL_HOME}/${HAIL_JAR}
#EOT
#
## Prepare bashrc and redownload script
#cat <<EOT >> /etc/bash.bashrc
#export HAIL_VERSION=${HAIL_VERSION}
#export SPARK_VERSION=${SPARK_VERSION}
#export SPARK_HOME=${SPARK_HOME}
#export HAIL_HOME=${HAIL_HOME}
#export HAIL_HASH=${HAIL_HASH}
#export HAIL_JAR=${HAIL_JAR}
#export HAIL_PYTHON_ZIP=${HAIL_PYTHON_ZIP}
#export _JAVA_OPTIONS='-Xmx8096m'
#export PYTHONPATH=/home/hail:${SPARK_HOME}/python:$(ls ${SPARK_HOME}/python/lib/py4j-*-src.zip):${HAIL_HOME}/${HAIL_PYTHON_ZIP}
#export SPARK_CONF_DIR=/home/hail/conf-ipython
#export PATH=/opt/conda/bin:$PATH
#EOT
#cat <<EOT >> /redownload_hail.sh
#mkdir -p $HAIL_HOME
#gsutil cp gs://hail-common/builds/${HAIL_VERSION}/jars/${HAIL_JAR} gs://hail-common/builds/${HAIL_VERSION}/python/${HAIL_PYTHON_ZIP} $HAIL_HOME
#EOT
#chmod +x /redownload_hail.sh
#
#cd /home/hail/gnomad_hail
#export PYTHONPATH=/home/hail:${SPARK_HOME}/python:$(ls ${SPARK_HOME}/python/lib/py4j-*-src.zip):${HAIL_HOME}/${HAIL_PYTHON_ZIP}
#export SPARK_CLASSPATH=${HAIL_HOME}/${HAIL_JAR}

# conda install -c damianavila82 rise
# conda install selenium phantomjs pillow -c conda-forge

#PYTHONPATH=.:$PYTHONPATH python tests/test_utils.py
# Use this instead to turn off Slack messages
# python -m unittest discover &> tests.log

# curl http://www.freefontspro.com/d/14454/arial.zip > /home/anaconda2/lib/python2.7/site-packages/matplotlib/mpl-data/fonts/ttf/arial.zip
# unzip /home/anaconda2/lib/python2.7/site-packages/matplotlib/mpl-data/fonts/ttf/arial.zip
