
"""
Module to expose more detailed version info for the installed `scipy`
"""
version = "1.17.0"
full_version = version
short_version = version.split('.dev')[0]
git_revision = "8c75ae75176236f233824e9a0483c26a69e6dfec"
release = 'dev' not in version and '+' not in version

if not release:
    version = full_version
