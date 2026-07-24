=======================
Vendored library policy
=======================

To simplify Bleach development, we're now vendoring certain libraries that
we use.

Vendored libraries must follow these rules:

1. Vendored libraries must be pure Python--no compiling.
2. Source code for the libary is included in this directory.
3. License must be included in this repo and in the Bleach distribution.
4. Requirements of the library become requirements of Bleach.
5. No modifications to the library may be made.


Adding/Updating a vendored library
==================================

Way to vendor a library or update a version:

1. Update ``vendor.txt`` with the library, version, and hash. You can use
   `hashin <https://pypi.org/project/hashin/>`_.
2. Remove all old files and directories of the old version.
3. Run ``pip_install_vendor.sh`` and check everything it produced in including
   the ``.dist-info`` directory and contents.
4. Update the bleach minor version in the next release.


Reviewing a change involving a vendored library
===============================================

Way to verify a vendored library addition/update:

1. Pull down the branch.
2. Delete all the old files and directories of the old version.
3. Run ``pip_install_vendor.sh``.
4. Run ``git diff`` and verify there are no changes.


NB: the current ``vendor.txt`` was generated with pip 20.2.3, which might be necessary to reproduce the dist-info


Removing/Unvendoring a vendored library
=======================================

A vendored library might be removed for any of the following reasons:

* it violates the vendoring policy (e.g. an incompatible license
  change)
* a suitable replacement is found
* bleach has the resources to test and QA new bleach releases against
  multiple versions of the previously vendored library

To unvendor a library:

1. Remove the library and its hashes from ``vendor.txt``.
2. Remove library files and directories from this directory.
3. Run ``install_vendor.sh`` and check the previously vendored library including
   the ``.dist-info`` directory and contents is not installed.
4. Update the bleach minor version in the next release.
