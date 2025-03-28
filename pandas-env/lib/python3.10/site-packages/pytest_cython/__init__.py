try:
    from pkg_resources import get_distribution
    __version__ = get_distribution('pytest-cython').version
    del get_distribution
except Exception:
    import warnings
    warnings.warn('could not get pytest-cython version; pkg_resources '
                  'not available or package not installed')
    __version__ = '0.0.0'
