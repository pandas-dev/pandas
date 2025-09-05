.. _develop.policies:

********
Policies
********

.. _policies.version:

Version policy
~~~~~~~~~~~~~~

pandas uses a loose variant of semantic versioning (`SemVer`_) to govern
deprecations, API compatibility, and version numbering.

A pandas release number is made up of ``MAJOR.MINOR.PATCH``.

API breaking changes should only occur in **major** releases. These changes
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
deprecation removed in the next major release (2.0.0).

.. note::

   pandas will sometimes make *behavior changing* bug fixes, as part of
   minor or patch releases. Whether or not a change is a bug fix or an
   API-breaking change is a judgement call. We'll do our best, and we
   invite you to participate in development discussion on the issue
   tracker or mailing list.

These policies do not apply to features marked as **experimental** in the documentation.
pandas may change the behavior of experimental features at any time.

.. _policies.python_support:

Python support
~~~~~~~~~~~~~~

pandas mirrors the `SPEC 0 guideline for Python support <https://scientific-python.org/specs/spec-0000>`__.

Security policy
~~~~~~~~~~~~~~~

To report a security vulnerability to pandas, please go to https://github.com/pandas-dev/pandas/security/policy and see the instructions there.

.. _SemVer: https://semver.org
