
"""
Module to expose more detailed version info for the installed `numpy`
"""
version = "2.4.1"
__version__ = version
full_version = version

git_revision = "d24bb7f48d3b0e3471c68f1309c130d0b65ee72a"
release = 'dev' not in version and '+' not in version
short_version = version.split("+")[0]
