"""Pytest configuration file for Hail tests."""

import pytest
import hail as hl


@pytest.fixture(scope="session", autouse=True)
def setup_hail():
    """
    Set up Hail before running tests.

    This fixture will ensure that the Hail context is initialized only once per test
    session and will be stopped after all tests are completed.
    """
    if not hl.utils.java.Env.is_fully_initialized():
        hl.init(log="/dev/null")

    yield
    hl.stop()
