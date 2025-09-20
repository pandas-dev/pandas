
"""
Module to expose more detailed version info for the installed `scipy`
"""
version = "1.16.2"
full_version = version
short_version = version.split('.dev')[0]
git_revision = "b1296b9b4393e251511fe8fdd3e58c22a1124899"
release = 'dev' not in version and '+' not in version

if not release:
    version = full_version
