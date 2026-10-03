# -*- coding: utf-8 -*-

"""
mbe_step
A SEAMM plug-in for many-body expansion (MBE) corrections of periodic cells and clusters
"""

# Bring up the classes so that they appear to be directly in
# the mbe_step package.

from .mbe import Mbe  # noqa: F401, E501
from .mbe_parameters import MbeParameters  # noqa: F401, E501
from .mbe_step import MbeStep  # noqa: F401, E501
from .tk_mbe import TkMbe  # noqa: F401, E501

from .metadata import metadata  # noqa: F401

# Handle versioneer
from ._version import get_versions

__author__ = "Paul Saxe"
__email__ = "psaxe@molssi.org"
versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
