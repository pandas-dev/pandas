# #############################################################################
# Copyright 2018 Hoffmann-La Roche
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# #############################################################################


from .pyreadstat import read_sav, read_sas7bdat, read_xport, read_dta, read_por, read_sas7bcat
from .pyreadstat import write_sav, write_dta, write_xport, write_por
from .pyreadstat import read_file_in_chunks, read_file_multiprocessing
from .pyclasses import metadata_container
from ._readstat_parser import ReadstatError, PyreadstatError
from .pyfunctions import set_value_labels, set_catalog_to_sas

__version__ = "1.3.5"

__all__ = (
    "read_sav",
    "read_sas7bdat",
    "read_xport",
    "read_dta",
    "read_por",
    "read_sas7bcat",
    "write_sav",
    "write_dta",
    "write_xport",
    "write_por",
    "read_file_in_chunks",
    "read_file_multiprocessing",
    "metadata_container",
    "ReadstatError",
    "PyreadstatError",
    "set_value_labels",
    "set_catalog_to_sas",
)
