"""
Shared pytest setup for the unit tests (tests/unit/test_<stage>.py, one file per stage).

The unit tests are fast (small images, seconds per file) and complement the stage scripts in
tests/, which run the slow brute-force checks and save figures. Run them all with

    python -m pytest tests/unit -q

Fixtures:
    in_tmp: run the test inside its own temporary folder (files the package writes, such as
        PSF caches, figures and parameter files, land there and are removed afterwards).

The tests use tiny images at 4 pixels per lambda/D, below the vortex model's MIN_SAMPLING, so its
SamplingWarning is silenced here (test_scene_stage6_outputs.py checks the warning itself).
"""
import os
import sys
import matplotlib

matplotlib.use("Agg")  # no windows; plt.show() becomes a no-op
import pytest

# tests/ holds helper modules shared with the stage scripts (extra_backends, extra_sources)
TESTS_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if TESTS_DIR not in sys.path:
    sys.path.insert(0, TESTS_DIR)


def pytest_configure(config):
    """Silence the vortex model's SamplingWarning (see the module docstring)."""
    config.addinivalue_line("filterwarnings", "ignore::toycoronagraph.optics.SamplingWarning")


@pytest.fixture
def in_tmp(tmp_path, monkeypatch):
    """Change the working directory to a fresh temporary folder for one test.

    Returns:
        pathlib.Path: The folder (also the working directory during the test).
    """
    monkeypatch.chdir(tmp_path)
    return tmp_path
