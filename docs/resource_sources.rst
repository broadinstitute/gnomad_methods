Resource Sources
================

gnomAD data is available through `multiple cloud providers' public datasets programs <https://gnomad.broadinstitute.org/news/2020-10-open-access-to-gnomad-data-on-multiple-cloud-providers/>`_.

The functions in the :doc:`gnomad.resources </api_reference/resources/index>` package can be configured to load data from different sources.

By default, resources are loaded from the gnomAD project's public Google Cloud Storage bucket ``gs://gnomad-public-requester-pays``. This bucket is `requester pays <https://cloud.google.com/storage/docs/requester-pays>`_, meaning that charges for data access and transfer will be billed to your Google Cloud project.

To load resources from a different source (for example, Google Public Datasets), use:

.. code-block:: python

    from gnomad.resources.config import gnomad_public_resource_configuration, GnomadPublicResourceSource

    gnomad_public_resource_configuration.source = GnomadPublicResourceSource.GOOGLE_CLOUD_PUBLIC_DATASETS

To see all available public sources for gnomAD resources, use:

.. code-block:: python

    from gnomad.resources.config import GnomadPublicResourceSource

    list(GnomadPublicResourceSource)


Custom Sources
--------------

Alternatively, instead of using one of the pre-defined public sources, a custom source can be provided.

.. code-block:: python

    from gnomad.resources.config import gnomad_public_resource_configuration

    gnomad_public_resource_configuration.source = "gs://my-bucket/gnomad-resources"
