# PDEP-17: Backwards compatibility and deprecation policy

- Created: 27 June 2024
- Status: Accepted
- Discussion: [#59125](https://github.com/pandas-dev/pandas/issues/59125)
- Author: [Abdulaziz Aloqeely](https://github.com/Aloqeely)
- Revision: 1

## Abstract

This PDEP defines pandas' backwards compatibility and deprecation policy.

The main additions to [pandas' current version policy](https://pandas.pydata.org/pandas-docs/version/2.2/development/policies.html) are:
- Deprecated functionality should remain unchanged in at least 2 minor releases before being changed or removed.
- Deprecations should initially use DeprecationWarning, and then be switched to FutureWarning in the last minor release before the major release they are planned to be removed in

## Motivation

Having a clear backwards compatibility and deprecation policy is crucial to having a healthy ecosystem. We want to ensure users can rely on pandas being stable while still allowing the library to evolve.

This policy will ensure that users have enough time to deal with deprecations while also minimizing disruptions on downstream packages' users.

## Scope

This PDEP covers pandas' approach to backwards compatibility and the deprecation and removal process.

## Background

pandas uses a loose variant of semantic versioning.  
A pandas release number is written in the format of ``MAJOR.MINOR.PATCH``.

## General policy

This policy applies to the [public API][1]. Anything not part of the [public API][1] or is marked as "Experimental" may be changed or removed at anytime.

- Breaking backwards compatibility should benefit more than it harms users.
- Breaking changes should go through a deprecation cycle before being implemented if possible.
- Breaking changes should only occur in major releases.
- No deprecations should be introduced in patch releases.
- Deprecated functionality should remain unchanged in at least 2 minor releases before being changed or removed.

Some bug fixes may require breaking backwards compatibility. In these cases, a deprecation cycle is not necessary. However, bug fixes which have a large impact on users might be treated as a breaking change. Whether or not a change is a bug fix or an API breaking change is a judgement call.

## Deprecation process

Deprecation provides a way to warn developers and give them time to adapt their code to the new functionality before the old behavior is eventually removed.

A deprecation's warning message should:
- Provide information on what is changing.
- Mention how to achieve similar behavior if an alternative is available.
- For large-scale deprecations, it is recommended to include a reason for the deprecation, alongside a discussion link to get user feedback.

Additionally, when one introduces a deprecation, they should:
- Use the appropriate warning class. More info on this can be found below.
- Add the GitHub issue/PR number as a comment above the warning line.
- Add an entry in the release notes.
- Mention that the functionality is deprecated in the documentation using the ``.. deprecated::`` directive.

### Which warning class to use

Starting in pandas 3.0, pandas will provide version-specific warnings. For example, ``Pandas4Warnings`` for all deprecation warnings that will be enforced in pandas 4.0. In addition to these, pandas exposes four additional classes to give users more control over pandas deprecation warnings.

- :class:`pandas.errors.PandasChangeWarning`: Base class of all pandas deprecation warnings.
- :class:`pandas.errors.PandasPendingDeprecationWarning`: Base class of all warnings which emit a ``PendingDeprecationWarning``, independent of the version they will be enforced.
- :class:`pandas.errors.PandasDeprecationWarning`: Base class of all warnings which emit a ``DeprecationWarning``, independent of the version they will be enforced.
- :class:`pandas.errors.PandasFutureWarning`: Base class of all warnings which emit a ``FutureWarning``, independent of the version they will be enforced.

### Enforcement of deprecations

When one enforces a deprecation, they should:
- Add an entry in the release notes.
- For API changes, replace the ``.. deprecated::`` directive in the documentation with a ``.. versionchanged::`` directive.

### PDEP-17 History

- 27 June 2024: Initial version.

[1]: https://pandas.pydata.org/docs/reference/index.html
