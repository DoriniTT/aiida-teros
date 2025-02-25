#!/usr/bin/env python

"""
Tests for the TEROS workflow.
"""

import os
import pytest
import yaml
from pathlib import Path

# Import functions to test
from aiida_teros.utils.io import load_yaml_config
from aiida_teros.schemas.config_schema import validate_config

# Test data directory
TEST_DATA_DIR = Path(__file__).parent / "test_data"


def test_config_validation():
    """Test that a valid config file passes validation."""
    config_path = TEST_DATA_DIR / "config.yaml"
    config = load_yaml_config(config_path)
    assert validate_config(config) is True


# Add more tests as needed