"""Pytest configuration file for Hail tests."""

import os

# CRITICAL: Force the local Spark backend BEFORE anything imports Hail.
#
# Hail reads `~/.config/hail/config.ini` at `hl.init()` time, which on
# developer machines is often configured to use the Hail Batch service
# backend (`[query] backend = batch`). When pytest collects test modules
# they `import hail as hl`, which can trigger Hail's lazy initialization
# during module collection — *before* any pytest fixture has had a chance
# to run. Hail then picks up the Batch backend from the global config and
# every subsequent Hail call submits jobs to Hail Batch, even though the
# `setup_hail` fixture below passes `backend="spark"` to `hl.init()`
# (which becomes a no-op because Hail is already initialized).
#
# Setting HAIL_QUERY_BACKEND=spark in the environment overrides
# `~/.config/hail/config.ini` at the source. Verified empirically: removing
# this line causes the test suite to start submitting Hail Batch jobs even
# with the explicit `backend="spark"` kwarg in the fixture.
os.environ.setdefault("HAIL_QUERY_BACKEND", "spark")

# Ensure Spark binds to a local loopback address before Hail/Spark are imported.
# On macOS the machine's hostname often does not resolve to a routable IP, which
# causes `sparkDriver` startup to fail with `BindException`. Setting
# `spark.driver.bindAddress` and `spark.driver.host` to 127.0.0.1 sidesteps the
# hostname resolution entirely. This must run before `import hail` because Hail
# constructs its SparkContext at init time using whatever PYSPARK_SUBMIT_ARGS
# was at process start. Verified empirically: removing this block causes
# `sparkDriver` startup to fail with a Py4JJavaError BindException.
_SPARK_BIND_ARGS = (
    "--conf spark.driver.bindAddress=127.0.0.1 " "--conf spark.driver.host=127.0.0.1"
)
_existing_args = os.environ.get("PYSPARK_SUBMIT_ARGS", "")
if "spark.driver.bindAddress" not in _existing_args:
    if _existing_args:
        os.environ["PYSPARK_SUBMIT_ARGS"] = f"{_SPARK_BIND_ARGS} {_existing_args}"
    else:
        os.environ["PYSPARK_SUBMIT_ARGS"] = f"{_SPARK_BIND_ARGS} pyspark-shell"

import hail as hl  # noqa: E402
import pytest  # noqa: E402


@pytest.fixture(scope="session", autouse=True)
def setup_hail():
    """
    Set up Hail before running tests.

    This fixture will ensure that the Hail context is initialized only once per
    test session and will be stopped after all tests are completed.

    Tests are forced to run on the local Spark backend so they never submit
    jobs to Hail Batch, regardless of the developer's global `hailctl`
    configuration. The actual backend selection is done at the top of this
    module via `HAIL_QUERY_BACKEND=spark` (which Hail reads at init time and
    which overrides `~/.config/hail/config.ini`). The `PYSPARK_SUBMIT_ARGS`
    setup at the top of this module pins Spark's driver to 127.0.0.1 to avoid
    hostname resolution issues on macOS.
    """
    if not hl.utils.java.Env.is_fully_initialized():
        hl.init(log="/dev/null")

    yield
    hl.stop()
