# Hail utilities for gnomAD

This repo contains a number of [Hail](https://hail.is/) utility functions and scripts for the [gnomAD project](http://gnomad.broadinstitute.org) and the [MacArthur lab](http://macarthurlab.org). As we continue to expand the size of our datasets, we are constantly seeking to find ways to reduce the complexity of our workflows and to make these functions more generic. As a result, the interface for many of these functions will change over time as we generalize their implementation for more flexible use within our scripts. We are also continously adapting our code to regular changes in the Hail interface. These repos thus represent only a snapshot of the gnomAD code base and are shared without guarantees or warranties.

We therefore encourage users to browse through the code and identify modules and functions that will be useful in their own pipelines, and to edit and reconfigure relevant code to suit their particular analysis and QC needs. 


Load using:

```python
from gnomad_hail.resources import *
from gnomad_hail.utils import *
from gnomad_hail.slack_utils import *
```

### resources.py

Resource map

### utils.py

Hail helper functions

### slack_utils.py

Helper functions for posting to Slack.
Submit along with `slack_creds.py` (e.g. in `HAIL_SCRIPTS` environment variable in `pyhail.py`) which has only a `slack_token` variable with a Slack API key.

### pyhail.py

Submission script for Hail. Notable differences from [cloudtools](http://github.com/nealelab/cloud-tools) include:

* Inline scripts (`--inline "print hc.read('gs://bucket/my.vds').count()"`)
* Add helper scripts designated in `HAIL_SCRIPTS` environment variable

### init_scripts

Google Cloud initialization scripts.
Cluster should be started up with `--init gs://gnomad-public/tools/inits/master-init.sh`, which is tracked here, deployed to cloud for usage.
`gnomad-init.sh` and `sparklyr-init.sh` are pulled from this repo and called by `master-init.sh`.

### tests

Tests are run automatically on creation of each cluster (in `gnomad-init.sh`) and can be triggered manually using: 
```bash
cd /path/to/gnomad_hail
PYTHONPATH=.:$PYTHONPATH python -m unittest discover
```
