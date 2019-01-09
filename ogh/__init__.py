from __future__ import (absolute_import,
                        division,
                        print_function,
                        unicode_literals)

import warnings
warnings.filterwarnings('ignore')
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
from .ogh import *
from . import ogh_xarray_landlab as oxl
from .ogh_meta import meta_file

__author__ = 'Jimmy Phuong'


class ogh_meta:
    """
    The json object that describes the Gridded climate data products
    """
    def __init__(self):
        self.__meta_data = dict(meta_file())

    # key-value retrieval
    def __getitem__(self, key):
        return(self.__meta_data[key])

    # key list
    def keys(self):
        return(self.__meta_data.keys())

    # value list
    def values(self):
        return(self.__meta_data.values())
