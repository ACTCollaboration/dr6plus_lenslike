"""
DR6+ Lensing Likelihood

Main likelihood class and helper functions for ACT DR6+ CMB lensing data.
"""

import numpy as np
import warnings
from typing import Dict, Optional, Tuple, Any
from scipy.interpolate import interp1d
import os

try:
    from cobaya.likelihoods.base_classes import InstallableLikelihood
except ImportError:
    # Fallback if cobaya is not installed
    InstallableLikelihood = object

#