
"""
Module to expose more detailed version info for the installed `numpy`
"""
version = "2.4.6"
__version__ = version
full_version = version

git_revision = "b832a09cf2a169c833dd2371e7c07aa00b293242"
release = 'dev' not in version and '+' not in version
short_version = version.split("+")[0]
