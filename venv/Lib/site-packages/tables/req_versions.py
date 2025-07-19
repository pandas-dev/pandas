"""Required versions for PyTables dependencies."""

from packaging.version import Version

# **********************************************************************
#   Runtime requirements, keep versions in sync with
#   "Prerequisites" in the User's Guide.
# **********************************************************************

# Minimum recommended versions for mandatory packages
min_numpy_version = Version("1.20.0")
min_numexpr_version = Version("2.6.2")
# These are library versions, not the python modules
min_hdf5_version = Version("1.10.5")
min_blosc_version = Version("1.11.1")
min_blosc2_version = Version("2.11.0")
