from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

# FIXME: Remove me once typed pointers are no longer supported.
def _opaque_pointers_enabled():
  import os
  return os.environ.get('LLVMLITE_ENABLE_OPAQUE_POINTERS', '0') == '1'
opaque_pointers_enabled = _opaque_pointers_enabled()
del _opaque_pointers_enabled
