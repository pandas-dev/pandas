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
will preserve the existing behavior while emitting a warning that provides
guidance on:

* How to achieve similar behavior if an alternative is available.
* The pandas version in which the deprecation will be enforced.

We will not introduce new deprecations in patch releases.

Deprecations will only be enforced in **major** releases. For example, if a
behavior is deprecated in pandas 3.1.0, it will continue to work, with a
warning, for all releases in the 3.x series. The behavior will change and the
deprecation removed in the next major release (4.0.0).

Deprecation warnings
^^^^^^^^^^^^^^^^^^^^

Following `PDEP-17`_, pandas uses a three-stage deprecation policy:

1. When a deprecation is first introduced in a minor release, it emits a
   :class:`~pandas.errors.PandasDeprecationWarning`. By default, Python hides
   ``DeprecationWarning`` outside the ``__main__`` module, so end users
   typically will not see these warnings unless they explicitly opt in. This
   gives downstream libraries time to migrate before their users are notified.
2. In the final minor release before the next major version, the warning is
   promoted to :class:`~pandas.errors.PandasFutureWarning`, which Python shows
   by default. End users see these warnings and have one minor release in
   which to update their code.
3. The deprecated behavior is removed in the next major release.

All warnings about upcoming changes inherit from
:class:`~pandas.errors.PandasChangeWarning`. Users can filter or surface
warnings using the following subclasses:

* :class:`~pandas.errors.Pandas4Warning` — changes that will be enforced in pandas 4.0.
* :class:`~pandas.errors.Pandas5Warning` — changes that will be enforced in pandas 5.0.
* :class:`~pandas.errors.PandasPendingDeprecationWarning` — emits a ``PendingDeprecationWarning``, regardless of target version.
* :class:`~pandas.errors.PandasDeprecationWarning` — emits a ``DeprecationWarning``, regardless of target version.
* :class:`~pandas.errors.PandasFutureWarning` — emits a ``FutureWarning``, regardless of target version.

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
.. _PDEP-17: https://pandas.pydata.org/pdeps/0017-backwards-compatibility-and-deprecation-policy.html
