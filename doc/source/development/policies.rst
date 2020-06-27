.. _develop.policies:

********
Policies
********

.. _policies.version:

Version policy
~~~~~~~~~~~~~~

.. versionchanged:: 1.0.0

pandas uses a loose variant of semantic versioning (`SemVer`_) to govern
deprecations, API compatibility, and version numbering.

A pandas release number is made up of ``MAJOR.MINOR.PATCH``.

API breaking changes should only occur in **major** releases. Theses changes
will be documented, with clear guidance on what is changing, why it's changing,
and how to migrate existing code to the new behavior.

Whenever possible, a deprecation path will be provided rather than an outright
breaking change.

pandas will introduce deprecations in **minor** releases. These deprecations
will preserve the existing behavior while emitting a warning that provide
guidance on:

* How to achieve similar behavior if an alternative is available
* The pandas version in which the deprecation will be enforced.

We will not introduce new deprecations in patch releases.

Deprecations will only be enforced in **major** releases. For example, if a
behavior is deprecated in pandas 1.2.0, it will continue to work, with a
warning, for all releases in the 1.x series. The behavior will change and the
deprecation removed in the next next major release (2.0.0).

.. note::

   pandas will sometimes make *behavior changing* bug fixes, as part of
   minor or patch releases. Whether or not a change is a bug fix or an
   API-breaking change is a judgement call. We'll do our best, and we
   invite you to participate in development discussion on the issue
   tracker or mailing list.

These policies do not apply to features marked as **experimental** in the documentation.
pandas may change the behavior of experimental features at any time.

Python support
~~~~~~~~~~~~~~

pandas will only drop support for specific Python versions (e.g. 3.6.x, 3.7.x) in
pandas **major** releases.

.. _SemVer: https://semver.org
