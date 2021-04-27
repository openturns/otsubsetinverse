"""
    otsubsetinverse --- An OpenTURNS module
    ==================================

    Contents
    --------
      'otsubsetinverse' is a module for OpenTURNS

"""

import sys
if sys.platform.startswith('win'):
    # this ensures OT dll is loaded
    import openturns

from .otsubsetinverse import *

__version__ = '1.8'

