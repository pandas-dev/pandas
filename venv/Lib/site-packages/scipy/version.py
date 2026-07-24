
"""
Module to expose more detailed version info for the installed `scipy`
"""
version = "1.18.0"
full_version = version
short_version = version.split('.dev')[0]
git_revision = "54ef5423f2e4376230ec3bfda6912a07a50958e3"
release = 'dev' not in version and '+' not in version

if not release:
    version = full_version
