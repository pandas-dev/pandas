try:
    import numexpr as ne
    _NUMEXPR_INSTALLED = True
    import distutils.version
    if distutils.version.LooseVersion(ne.__version__) < '2':
        _NUMEXPR_INSTALLED = False
        import warnings
        warnings.warn("pandas requires at least numexpr version 2.0"
                      " to accelerate with numexpr. Found numexpr version"
                      " %s" % ne.__version__)
    del ne
    del distutils.version
except ImportError:  # pragma: no cover
    _NUMEXPR_INSTALLED = False

#: tracks globally whether numexpr should be used.
_USE_NUMEXPR = _NUMEXPR_INSTALLED
