""" support pyarrow compatibility across versions """

from distutils.version import LooseVersion

try:
    import pyarrow as pa

    _pa_version = pa.__version__
    _palv = LooseVersion(_pa_version)
    pa_version_under1p0 = _palv < LooseVersion("1.0.0")
    pa_version_under2p0 = _palv < LooseVersion("2.0.0")
    pa_version_under3p0 = _palv < LooseVersion("3.0.0")
    pa_version_under4p0 = _palv < LooseVersion("4.0.0")
except ImportError:
    pa_version_under1p0 = True
    pa_version_under2p0 = True
    pa_version_under3p0 = True
    pa_version_under4p0 = True
