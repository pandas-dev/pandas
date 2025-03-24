########################################################################
#
#       Created: September 22, 2010
#       Author:  The Blosc development team - blosc@blosc.org
#
########################################################################


from blosc.version import __version__

# Blosc C symbols that we want to export
from blosc.blosc_extension import (
    BLOSC_VERSION_STRING as VERSION_STRING,
    BLOSC_VERSION_DATE as VERSION_DATE,
    BLOSC_MAX_BUFFERSIZE as MAX_BUFFERSIZE,
    BLOSC_MAX_THREADS as MAX_THREADS,
    BLOSC_MAX_TYPESIZE as MAX_TYPESIZE,
    BLOSC_NOSHUFFLE as NOSHUFFLE,
    BLOSC_SHUFFLE as SHUFFLE,
    BLOSC_BITSHUFFLE as BITSHUFFLE,
    init,
    destroy,
    )

# Restored old symbols for backward compatibility with pre 1.3.0
BLOSC_VERSION_STRING = VERSION_STRING
BLOSC_VERSION_DATE = VERSION_DATE
BLOSC_MAX_BUFFERSIZE = MAX_BUFFERSIZE
BLOSC_MAX_THREADS = MAX_THREADS
BLOSC_MAX_TYPESIZE = MAX_TYPESIZE

# Translation of filters to strings
filters = {NOSHUFFLE: "noshuffle",
           SHUFFLE: "shuffle",
           BITSHUFFLE: "bitshuffle"}

from blosc.toplevel import (
    compress,
    compress_ptr,
    decompress,
    decompress_ptr,
    pack_array,
    unpack_array,
    detect_number_of_cores,
    free_resources,
    set_nthreads,
    set_blocksize,
    get_blocksize,
    set_releasegil,
    compressor_list,
    code_to_name,
    name_to_code,
    clib_info,
    get_clib,
    get_cbuffer_sizes,
    cbuffer_validate,
    print_versions,
    )

# Dictionaries for the maps between compressor names and libs
cnames = compressor_list()
# Map for compression names and libs
cname2clib = {name: clib_info(name)[0] for name in cnames}
# Map for compression libraries and versions
clib_versions = dict(clib_info(name) for name in cnames)


# Initialize Blosc
init()
# default to keep GIL, since it's just extra overhead if we aren't
# threading ourselves
set_releasegil(False)
# Internal Blosc threading
nthreads = ncores = detect_number_of_cores()
# Protection against too many cores
nthreads = min(nthreads, 8)
set_nthreads(nthreads)
blosclib_version = "%s (%s)" % (VERSION_STRING, VERSION_DATE)
import atexit
atexit.register(destroy)


# Tests
from blosc.test import run as test
