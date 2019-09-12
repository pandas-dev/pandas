.. _develop.policies:

********
Policies
********

.. _policies.version:

Version Policy
~~~~~~~~~~~~~~

.. versionchanged:: 1.0.0

Pandas uses a version of `SemVer`_ to govern deprecations, API compatibility, and version numbering.

A pandas release number is made up of ``MAJOR.MINOR.PATCH``.

API breaking changes should only occur in **major** releases. Theses changes will be documented,
with clear guidance on what is changing, why it's changing, and how to migrate existing code to the
new beahvior.

Whenever possible, a deprecation path will be provided rather than an outright breaking change.

Pandas will introduce deprecations in **major** or **minor** releases. These deprecations will
preserve the existing behavior while emitting a warning that provide guidance

* How to achieve similar behavior if an alternative is available
* The pandas version in which the deprecation will be enforced.

We will not introduce new deprecations in patch releases.

Deprecations will only be enforced in **major** releases.

.. note::

   Pandas will sometimes make *behavior changing* bug fixes, as part of
   minor or patch releases. Whether or not a change is a bug fix or an
   API-breaking change is a judgement call. We'll do our best, and we
   invite you to participate in development discussion on the issue
   tracker or mailing list.


.. _SemVer: https://semver.org
