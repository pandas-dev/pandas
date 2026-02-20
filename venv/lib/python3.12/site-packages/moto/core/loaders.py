"""Module for loading various model files along with any `moto-extras`

Heavily leverages the existing botocore loading mechanism, which expects
a particular directory layout and file naming convention:

    <MOTO_ROOT>
      |
      |-- service_name
      |   |-- YYYY-MM-DD (matching API version for service_name)
      |       |-- service-2.moto-extras.json
      |       |-- paginators-1.moto-extras.json

The ``moto-extras`` files represent extra data to be applied to the model
after it is loaded. Data in these files could augment or modify the AWS
service model by supplying additional error shapes or default input
values for use by Moto.
"""

import os

from botocore.loaders import Loader as BotocoreLoader

import moto

moto_root = os.path.dirname(os.path.abspath(moto.__file__))


class Loader(BotocoreLoader):
    BUILTIN_EXTRAS_TYPES = ["moto"]


def create_loader() -> Loader:
    return Loader(extra_search_paths=[moto_root])
