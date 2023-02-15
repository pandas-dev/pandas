# PDEP-1: Purpose and guidelines

- Created: 3 August 2022
- Status: Accepted
- Discussion: [#47444](https://github.com/pandas-dev/pandas/pull/47444)
- Author: [Marc Garcia](https://github.com/datapythonista)
- Revision: 1

## PDEP definition, purpose and scope

A PDEP (pandas enhancement proposal) is a proposal for a **major** change in
pandas, in a similar way as a Python [PEP](https://peps.python.org/pep-0001/)
or a NumPy [NEP](https://numpy.org/neps/nep-0000.html).

Bug fixes and conceptually minor changes (e.g. adding a parameter to a function) are out of the
scope of PDEPs. A PDEP should be used for changes that are not immediate and not obvious, when
everybody in the pandas community needs to be aware of the possibility of an upcoming change.
Such changes require detailed documentation before being implemented and frequently lead to a
significant discussion within the community.

PDEP are appropriate for user facing changes, internal changes and significant discussions.
Examples of topics worth a PDEP could include substantial API changes, breaking behavior changes,
moving a module from pandas to a separate repository, a refactoring of the pandas block manager
or a proposal of a new code of conduct. It is not always trivial to know which issue has enough
scope to require the full PDEP process. Some simple API changes have sufficient consensus among
the core team, and minimal impact on the community. On the other hand, if an issue becomes
controversial, i.e. it generated a significant discussion, one could suggest opening a PDEP to
formalize and document the discussion, making it easier for the wider community to participate.
For context, see the list of issues that could have been a PDEP at the bottom of this page.

## PDEP guidelines

### Target audience

A PDEP is a public document available to anyone, but the main stakeholders to
consider when writing a PDEP are:

- The core development team, who will have the final decision on whether a PDEP
  is approved or not
- Contributors to pandas and other related projects, and experienced users. Their
  feedback is highly encouraged and appreciated, to make sure all points of views
  are taken into consideration
- The wider pandas community, in particular users, who may or may not have feedback
  on the proposal, but should know and be able to understand the future direction of
  the project

### PDEP authors

Anyone can propose a PDEP, but core members need to sponsor a proposal made by non-core
contributors. To submit a PDEP as a community member, please propose the PDEP concept on
[an issue](https://github.com/pandas-dev/pandas/issues/new/choose), and find a pandas team
member to collaborate with. They can co-author the PDEP with you and submit it to the PDEPs
repository.

### Workflow

The possible states of a PDEP are:

- Under discussion
- Accepted
- Implemented
- Rejected

Next is described the workflow that PDEPs can follow.

#### Submitting a PDEP

Proposing a PDEP is done by creating a PR adding a new file to `web/pdeps/`.
The file is a markdown file, you can use `web/pdeps/0001.md` as a reference
for the expected format.

The initial status of a PDEP will be `Status: Draft`. Once it is ready for discussion, the author(s)
change it to `Status: Under discussion`, and the following are notified: core and triage teams
and the pandas-dev mailing list. This will be changed to `Status: Accepted` when the PDEP is ready
and have the approval of the core team.


#### Accepted PDEP

A PDEP can only be accepted by the core development team, if the proposal is considered
worth implementing. Decisions will be made based on the process detailed in the
[pandas governance document](https://github.com/pandas-dev/pandas-governance/blob/master/governance.md).
In general, more than one approval will be needed before the PR is merged. And
there should not be any `Request changes` review at the time of merging.

Once a PDEP is accepted, any contributions can be made toward the implementation of the PDEP,
with an open-ended completion timeline. Development of pandas is difficult to understand and
forecast, being that the contributors to pandas are a mix of volunteers and developers paid from different sources,
with different priorities. For companies, institutions or individuals with interest in seeing a
PDEP being implemented, or to in general see progress to the pandas roadmap, please check how
you can help in the [contributing page](/contribute.html).

#### Implemented PDEP

Once a PDEP is implemented and available in the main branch of pandas, its
status will be changed to `Status: Implemented`, so there is visibility that the PDEP
is not part of the roadmap and future plans, but a change that has already
happened. The first pandas version in which the PDEP implementation is
available will also be included in the PDEP header with for example
`Implemented: v2.0.0`.

#### Rejected PDEP

A PDEP can be rejected when the final decision is that its implementation is
not in the best interests of the project. Rejected PDEPs are as useful as accepted
PDEPs, since there are discussions that are worth having, and decisions about
changes to pandas being made. They will be merged with `Status: Rejected`, so
there is visibility on what was discussed and what was the outcome of the
discussion. A PDEP can be rejected for different reasons, for example good ideas
that aren't backward-compatible, and the breaking changes aren't considered worth
implementing.

#### Invalid PDEP

For submitted PDEPs that do not contain proper documentation, are out of scope, or
are not useful to the community for any other reason, the PR will be closed after
discussion with the author, instead of merging them as rejected. This is to avoid
adding noise to the list of rejected PDEPs, which should contain documentation as
good as an accepted PDEP, but where the final decision was to not implement the changes.

## Evolution of PDEPs

Most PDEPs aren't expected to change after accepted. Once there is agreement in the changes,
and they are implemented, the PDEP will be only useful to understand why the development happened,
and the details of the discussion.

But in some cases, a PDEP can be updated. For example, a PDEP defining a procedure or
a policy, like this one (PDEP-1). Or cases when after attempting the implementation,
new knowledge is obtained that makes the original PDEP obsolete, and changes are
required. When there are specific changes to be made to the original PDEP, this will
be edited, its `Revision: X` label will be increased by one, and a note will be added
to the `PDEP-N history` section. This will let readers understand that the PDEP has
changed and avoid confusion.

## List of issues that could have been PDEPs for context
### Clear examples for potential PDEPs:
- Adding a new parameter to many existing methods, or deprecating one in many places. For example:
  - The numeric_only deprecation affected many methods and could have been a PDEP.
- Adding a new data type has impact on a variety of places that need to handle the data type.
  Such wide-ranging impact would require a PDEP. For example:
  - Categorical (GH-7217, GH-8074), StringDtype (GH-8640), ArrowDtype
- A significant (breaking) change in existing behavior. For example:
  - copy/view changes (GH-36195)
- Support of new python features with a wide impact on the project. For example:
  - Supporting typing within pandas vs. creation of pandas-stubs (GH-43197, GH-45253) New
    required dependency

### Borderline examples:
- Small changes to core functionality, such as DataFrame and Series, should always be considered
  as a PDEP candidate as it will likely have a big impact on users. But the same types of changes in other methods would not be good PDEP candidates. That said, any discussion, no matter how small the change, which becomes controversial is a PDEP candidate. Consider if more attention and/or a formal decision-making process would help. Following are some examples we hope can help clarify our meaning here:
- API breaking changes, or discussion thereof, could be a PDEP. For example:
  - Value counts rename (GH-49497). The scope doesn’t justify a PDEP at first, but later a
    discussion about whether it should be executed as breaking change or with deprecation emerges, which could benefit from a PDEP process.
- Adding new parameters or methods to an existing method typically won't require a PDEP for
  non-core features. For example:
  - Both dropna(percentage) (GH-35299), and Timestamp.normalize() (GH-8794) wouldn't have
    required a PDEP.
  - On the other hand, DataFrame.assign() might. While it's a single method without backwards
    compatibility concerns, it’s also a core feature and the discussion should be highly visible.
- Deprecating or removing a single method would not require a PDEP in most cases. For example:
  - DataFrame.xs (GH-6249) is an example of depredations on core features that would be good a
    candidate for a PDEP.
- Changing the default value of parameters in a core pandas method is another edge case. For
  example:
  - Such changes in dropna,  DataFrame.groupby, or in Series.groupby could be PDEPs.
- New top level modules and/or exposing internal classes. For example:
  - Add pandas.api.typing (GH-48577) is relatively small and wouldn’t necessarily require a PDEP.

### A counter-example:
- Significant changes to developer/contributor processes don't require a PDEP as they aren't
going to have an impact on users. For example:
  - Changing the build system to meson


### PDEP-1 History

- 3 August 2022: Initial version
- 15 February 2023: Revision 1
