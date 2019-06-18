"""
Various decorators to validate input/output arguments to functions.
"""
__all__ = []

import functools
import inspect
import numpy as np

from astropy import units as u
from plasmapy.utils.decorators.helpers import preserve_signature
from typing import (Any, Dict, List, Union)
