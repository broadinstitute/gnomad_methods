# Hail utilities for gnomAD

This repo contains a number of [Hail](https://hail.is/) utility functions and scripts for the [gnomAD project](http://gnomad.broadinstitute.org) and the [MacArthur lab](http://macarthurlab.org).
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
