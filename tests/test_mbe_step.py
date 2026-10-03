#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `mbe_step` package."""

import pytest  # noqa: F401
import mbe_step  # noqa: F401


def test_construction():
    """Just create an object and test its type."""
    result = mbe_step.Mbe()
    assert str(type(result)) == "<class 'mbe_step.mbe.Mbe'>"
