"""
This module contains models and functions to manage the library/list of models
"""
from __future__ import absolute_import
from . _libfunc import *

__all__ = ["updateModelList", "getModelNames", "getModelDesc"]

try:
    from . _modellist import *
    from . _modellist import _model_list

except ImportError:
    _model_list = None

if _model_list is not None:
    for i in range(len(_model_list)):
        __all__.append(_model_list[i])
