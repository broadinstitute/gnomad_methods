Getting Started
===============

This short guide will help you get started with hail and the gnomaD Python package. Please note this guide was created in 2022 and the arguments for creating hail clusters may be dated. For example, if you are using hail version 0.2.131 or later and dataproc version 2.2 or later, the command in step 2 below requires the flag `--public-ip-address` at the time this guide was updated, September 2024. Please refer to the `Hail documentation <https://hail.is/docs/0.2/getting_started.html>`_ for the most up-to-date information.


1. `Install Hail <https://hail.is/docs/0.2/getting_started.html#installation>`_::

    pip install hail

2. Use ``hailctl`` to start a `Google Dataproc <https://cloud.google.com/dataproc/>`_ cluster with the
   ``gnomad`` package installed (see `Hail on the Cloud <https://hail.is/docs/0.2/hail_on_the_cloud.html>`_ for more detail on ``hailctl``)::

    hailctl dataproc start cluster-name --packages gnomad

3. Connect to a `Jupyter Notebook <https://jupyter-notebook.readthedocs.io/en/stable/notebook.html>`_ on the cluster::

    hailctl dataproc connect cluster-name notebook

4. Import gnomAD data in `Hail Table <https://hail.is/docs/0.2/hail.Table.html>`_ format:

    * gnomAD v2.1.1 variants::

        from gnomad.resources.grch37 import gnomad

        gnomad_v2_exomes = gnomad.public_release("exomes")
        exomes_ht = gnomad_v2_exomes.ht()
        exomes_ht.describe()

        gnomad_v2_genomes = gnomad.public_release("genomes")
        genomes_ht = gnomad_v2_genomes.ht()
        genomes_ht.describe()

    * gnomAD v3 variants::

        from gnomad.resources.grch38 import gnomad
        gnomad_v3_genomes = gnomad.public_release("genomes")
        ht = gnomad_v3_genomes.ht()
        ht.describe()

5. Shut down the cluster when finished with it::

    hailctl dataproc stop cluster-name
