# PDEP-17: Backwards compatibility and deprecation policy

- Created: 27 June 2024
- Status: Under discussion
- Discussion: [#59125](https://github.com/pandas-dev/pandas/issues/59125)
- Author: [Abdulaziz Aloqeely](https://github.com/Aloqeely)
- Revision: 1

## Abstract

This PDEP defines pandas' backwards compatibility and deprecation policy.

## Motivation

Having a clear backwards compatibility and deprecation policy is crucial to having a healthy ecosystem. We want to ensure users can rely on pandas being stable while still allowing the library to evolve.

This policy will ensure that users have enough time to deal with deprecations while also minimizing disruptions on downstream packages' users.

## Scope

This PDEP covers pandas' approach to backwards compatibility and the deprecation and removal process.

The decision making process for deprecating features is out of scope for this PDEP.

## Background

pandas uses a loose variant of semantic versioning.  
A pandas release number is written in the format of ``MAJOR.MINOR.PATCH``.

For the purposes of this PDEP, the last minor release before a major release will be referred to as the 'final minor version'.

## General policy

This policy applies to the [public API][1]. Anything not part of the [public API][1] or is marked as "Experimental" may be changed or removed at anytime.

- Breaking backwards compatibility should benefit more than it harms users.
- Breaking changes should go through a deprecation cycle before being implemented.
- Breaking changes should only occur in major releases.
- No deprecations should be introduced in patch releases.

Some bug fixes may require breaking backwards compatibility. In these cases, a deprecation cycle is not necessary. However, bug fixes which have a large impact on users might be treated as a breaking change. Whether or not a change is a bug fix or an API breaking change is a judgement call.

## Deprecation process

Deprecation provides a way to warn developers and give them time to adapt their code to the new functionality before the old behavior is eventually removed.

A deprecation's warning message should:
- Provide information on what is changing.
- Mention how to achieve similar behavior if an alternative is available.
- Include the version in which the deprecation will be enforced.
- For large-scale deprecations, it is recommended to include a reason for the deprecation, alongside a discussion link to get user feedback.

Additionally, when one introduces a deprecation, they should:
- Use the appropriate warning class. More info on this can be found below.
- Add the GitHub issue/PR number as a comment above the warning line.
- Add an entry in the release notes.
- Mention that the functionality is deprecated in the documentation using the ``.. deprecated::`` directive.

### Which warning class to use

Deprecations should initially use ``DeprecationWarning``, and then be switched to ``FutureWarning`` for broader visibility in the final minor version before the major release they are planned to be removed in.  
This implementation detail can be ignored by using the appropriate ``PandasDeprecationWarning`` variable, which will be aliased to the proper warning class based on the pandas version.

Not all deprecations have to use ``DeprecationWarning`` but all deprecations should eventually transition to ``FutureWarning``, i.e. deprecations in the final minor version which are planned to be removed in the major release after will immediately use ``FutureWarning``.

It is recommended to not introduce large-scale deprecations in the final minor version which are planned to be removed in the major release after, since that will immediately be using a loud ``FutureWarning`` with not much time between deprecation and removal. Instead, a ``DeprecationWarning`` should be used, and the removal should be scheduled for a later major release.

### Support window for the final minor version

Defining a release policy is not part of this PDEP's scope, but, the final minor version plays a special role in deprecations, since:
- It is the version where deprecations planned for removal in the next major release get switched from a ``DeprecationWarning`` to a more visible ``FutureWarning``.
- ``FutureWarning`` deprecations released in it will immediately be removed in the next major release.

And so, this document recommends a minimum of 6 months between the final minor version and the major release after it.

### Enforcement of deprecations

When one enforces a deprecation, they should:
- Add an entry in the release notes.
- Replace the ``.. deprecated::`` directive in the documentation with a ``.. versioncchanged::`` directive.

### PDEP-17 History

- 27 June 2024: Initial version.

[1]: https://pandas.pydata.org/docs/reference/index.html
