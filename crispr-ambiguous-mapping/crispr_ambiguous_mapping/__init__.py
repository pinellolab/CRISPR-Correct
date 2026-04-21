#!/usr/bin/env python

# §7.6: library convention — attach a NullHandler to the package logger so
# that consumers who don't configure logging don't see a "No handlers could
# be found" warning. Users who want to see INFO output call
# `logging.basicConfig(level=logging.INFO)` once at startup.
import logging as _logging
_logging.getLogger(__name__).addHandler(_logging.NullHandler())
del _logging

from . import mapping
from . import utility
from . import models
from . import processing
from . import visualization
from . import quality_control
from . import postprocessing
