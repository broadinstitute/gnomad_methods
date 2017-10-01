# Hail utilities for gnomAD

This repo contains a number of Hail utility functions and scripts for the [gnomAD project](http://gnomad.broadinstitute.org) and the [MacArthur lab](http://macarthurlab.org).
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