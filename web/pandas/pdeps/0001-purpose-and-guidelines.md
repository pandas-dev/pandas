# PDEP-1: Purpose and guidelines

- Created: 3 August 2022
- Status: Accepted
- Discussion: [#47444](https://github.com/pandas-dev/pandas/pull/47444),
  [#51417](https://github.com/pandas-dev/pandas/pull/51417)
- Author: [Marc Garcia](https://github.com/datapythonista),
  [Noa Tamir](https://github.com/noatamir)
- Revision: 3

[TOC]

## PDEP definition, purpose and scope

A PDEP (pandas enhancement proposal) is a proposal for a **major** change in
pandas, in a similar way as a Python [PEP](https://peps.python.org/pep-0001/)
or a NumPy [NEP](https://numpy.org/neps/nep-0000.html).

Bug fixes and conceptually minor changes (e.g. adding a parameter to a function) are out of the
scope of PDEPs. A PDEP should be used for changes that are not immediate and not obvious, when
everybody in the pandas community needs to be aware of the possibility of an upcoming change.
Such changes require detailed documentation before being implemented and frequently lead to a
significant discussion within the community.

PDEPs are appropriate for user facing changes, internal changes and significant discussions.
Examples of topics worth a PDEP could include substantial API changes, breaking behavior changes,
moving a module from pandas to a separate repository, or a refactoring of the pandas block manager.
It is not always trivial to know which issue has enough scope to require the full PDEP process.
Some simple API changes have sufficient consensus among the core team, and minimal impact on the
community. On the other hand, if an issue becomes controversial, i.e. it generated a significant
discussion, one could suggest opening a PDEP to formalize and document the discussion, making it
easier for the wider community to participate. For context, see
[the list of issues that could have been a PDEP](#List-of-issues).

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

Anyone can propose a PDEP, but a core member should be engaged to advise on a proposal made by
non-core contributors. To submit a PDEP as a community member, please propose the PDEP concept on
[an issue](https://github.com/pandas-dev/pandas/issues/new/choose), and find a pandas team
member to collaborate with. They can advise you on the PDEP process and should be listed as an
advisor on the PDEP when it is submitted to the PDEP repository.

### Workflow

#### Rationale

Our workflow was created to support and enable a consensus seeking process, and to provide clarity,
for current and future authors, as well as voting members. It is not a strict policy, and we
discourage any interpretation which seeks to take advantage of it in a way that could "force" or
"sneak" decisions in one way or another. We expect and encourage transparency, active discussion,
feedback, and compromise from all our community members.

#### PDEP States

The possible states of a PDEP are:

- Draft
- Under discussion
- Accepted
- Implemented
- Rejected
- Withdrawn

Next is described the workflow that PDEPs can follow.

#### Submitting a PDEP

Proposing a PDEP is done by creating a PR adding a new file to `web/pandas/pdeps/`.
The file is a markdown file, you can use `web/pandas/pdeps/0001-purpose-and-guidelines.md` as a reference
for the expected format.

The initial status of a PDEP will be `Status: Draft`. This will be changed to
`Status: Under discussion` by the author(s), when they are ready to proceed with the decision
making process.

#### PDEP Discussion Timeline

A PDEP discussion will remain open for up to 60 days. This period aims to enable participation
from volunteers, who might not always be available to respond quickly, as well as provide ample
time to make changes based on suggestions and considerations offered by the participants.
Similarly, the following voting period will remain open for 15 days.

To enable and encourage discussions on PDEPs, we follow a notification schedule. At each of the
following steps, the pandas team, and the pandas-dev mailing list are notified via GitHub and
E-mail:

- Once a PDEP is ready for discussion.
- After 30 days, with a note that there is at most 30 days remaining for discussion,
  and that a vote will be called for if no discussion occurs in the next 15 days.
- After 45 days, with a note that there is at most 15 days remaining for discussion,
  and that a vote will be called for in 15 days.
- Once the voting period starts, after 60 days or in case of an earlier vote, with 15 days
  remaining for voting.
- After 10 voting days, with 5 days remaining for voting.

After 30 discussion days, in case 15 days passed without any new unaddressed comments,
the authors may close the discussion period preemptively, by sending an early reminder
of 15 days remaining until the voting period starts.

#### Casting Votes

As the voting period starts, a VOTE issue is created which links to the PDEP discussion pull request.
Each voting member, including author(s) with voting rights, may cast a vote by adding one of the following comments:

- +1: approve.
- 0: abstain.
  - Reason: A one sentence reason is required.
- -1: disapprove
  - Reason: A one sentence reason is required.

A disapprove vote requires prior participation in the PDEP discussion issue.

Comments made on the public VOTE issue by non-voting members will be deleted.

Once the voting period ends, any voter may tally the votes in a comment, using the format: w-x-y-z,
where w stands for the total of approving, x of abstaining, y of disapproving votes cast, and z
of number of voting members who did not respond to the VOTE issue. The tally of the votes will state
if a quorum has been reached or not.

#### Quorum and Majority

For a PDEP vote to result in accepting the proposal, a quorum is required. All votes (including
abstentions) are counted towards the quorum. The quorum is computed as the lower of these two
values:

- 11 voting members.
- 50% of voting members.

Given a quorum, a majority of 70% of the non-abstaining votes is required as well, i.e. 70% of
the approving and disapproving votes must be in favor.

Thus, abstaining votes count towards a quorum, but not towards a majority. A voting member might
choose to abstain when they have participated in the discussion, have some objections to the
proposal, but do not wish to stop the proposal from moving forward, nor indicate their full
support.

If a quorum was not reached by the end of the voting period, the PDEP is not accepted. Its status
will change to rejected.

#### Accepted PDEP

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
that are not backward-compatible, and the breaking changes are not considered worth
implementing.

The PDEP author(s) can also decide to withdraw the PDEP before a final decision
is made (`Status: Withdrawn`), when the PDEP authors themselves have decided
that the PDEP is actually a bad idea, or have accepted that it is not broadly
supported or a competing proposal is a better alternative.

The author(s) may choose to resubmit a rejected or withdrawn PDEP. We expect
authors to use their judgement in that case, as to whether they believe more
discussion, or an amended proposal has the potential to lead to a different
result. A new PDEP is then created, which includes a link to the previously
rejected PDEP.

#### Invalid PDEP

For submitted PDEPs that do not contain proper documentation, are out of scope, or
are not useful to the community for any other reason, the PR will be closed after
discussion with the author, instead of merging them as rejected. This is to avoid
adding noise to the list of rejected PDEPs, which should contain documentation as
good as an accepted PDEP, but where the final decision was to not implement the changes.

## Evolution of PDEPs

Most PDEPs are not expected to change after they are accepted. Once there is agreement on the changes,
and they are implemented, the PDEP will be only useful to understand why the development happened,
and the details of the discussion.

But in some cases, a PDEP can be updated. For example, a PDEP defining a procedure or
a policy, like this one (PDEP-1). Or cases when after attempting the implementation,
new knowledge is obtained that makes the original PDEP obsolete, and changes are
required. When there are specific changes to be made to the original PDEP, this will
be edited, its `Revision: X` label will be increased by one, and a note will be added
to the `PDEP-N history` section. This will let readers understand that the PDEP has
changed and avoid confusion.

## <a id="List-of-issues"></a> List of issues that could have been PDEPs for context
### Clear examples for potential PDEPs:

- Adding a new parameter to many existing methods, or deprecating one in many places. For example:
  - The `numeric_only` deprecation ([GH-28900][28900]) affected many methods and could have been a PDEP.
- Adding a new data type has impact on a variety of places that need to handle the data type.
  Such wide-ranging impact would require a PDEP. For example:
  - `Categorical` ([GH-7217][7217], [GH-8074][8074]), `StringDtype` ([GH-8640][8640]), `ArrowDtype`
- A significant (breaking) change in existing behavior. For example:
  - Copy/view changes ([GH-36195][36195])
- Support of new Python features with a wide impact on the project. For example:
  - Supporting typing within pandas vs. creation of `pandas-stubs` ([GH-43197][43197],
    [GH-45253][45253])
- New required dependency.
- Removing module from the project or splitting it off to a separate repository:
  - Moving rarely used I/O connectors to a separate repository [GH-28409](28409)
- Significant changes to contributors' processes are not going to have an impact on users, but
they do benefit from structured discussion among the contributors. For example:
  - Changing the build system to meson ([GH-49115][49115])

### Borderline examples:
Small changes to core functionality, such as `DataFrame` and `Series`, should always be
considered as a PDEP candidate as it will likely have a big impact on users. But the same types
of changes in other functionalities would not be good PDEP candidates. That said, any discussion,
no matter how small the change, which becomes controversial is a PDEP candidate. Consider if more
attention and/or a formal decision-making process would help. Following are some examples we
hope can help clarify our meaning here:

- API breaking changes, or discussion thereof, could be a PDEP. For example:
  - `value_counts` result rename ([GH-49497][49497]). The scope does not justify a PDEP at first, but later a
    discussion about whether it should be executed as a breaking change or with deprecation
    emerges, which could benefit from the PDEP process.
- Adding new methods or parameters to an existing method typically will not require a PDEP for
  non-core features. For example:
  - Both `dropna(percentage)` ([GH-35299][35299]), and `Timestamp.normalize()` ([GH-8794][8794])
    would not have required a PDEP.
  - On the other hand, `DataFrame.assign()` might. While it is a single method without backwards
    compatibility concerns, it is also a core feature and the discussion should be highly visible.
- Deprecating or removing a single method would not require a PDEP in most cases.
  - That said, `DataFrame.append` ([GH-35407][35407]) is an example of deprecations on core
    features that would be a good candidate for a PDEP.
- Changing the default value of parameters in a core pandas method is another edge case. For
  example:
  - Such changes for `dropna` in `DataFrame.groupby` and `Series.groupby` could be a PDEP.
- New top level modules and/or exposing internal classes. For example:
  - Add `pandas.api.typing` ([GH-48577][48577]) is relatively small and would not necessarily
    require a PDEP.


### PDEP-1 History

- 3 August 2022: Initial version ([GH-47938][47938])
- 15 February 2023: Version 2 ([GH-51417][51417]) clarifies the scope of PDEPs and adds examples
- 09 June 2023: Version 3 ([GH-53576][53576]) defines a structured decision making process for PDEPs

[7217]: https://github.com/pandas-dev/pandas/pull/7217
[8074]: https://github.com/pandas-dev/pandas/issues/8074
[8640]: https://github.com/pandas-dev/pandas/issues/8640
[36195]: https://github.com/pandas-dev/pandas/issues/36195
[43197]: https://github.com/pandas-dev/pandas/issues/43197
[45253]: https://github.com/pandas-dev/pandas/issues/45253
[49497]: https://github.com/pandas-dev/pandas/issues/49497
[35299]: https://github.com/pandas-dev/pandas/issues/35299
[8794]: https://github.com/pandas-dev/pandas/issues/8794
[6249]: https://github.com/pandas-dev/pandas/issues/6249
[48577]: https://github.com/pandas-dev/pandas/issues/48577
[49115]: https://github.com/pandas-dev/pandas/pull/49115
[28409]: https://github.com/pandas-dev/pandas/issues/28409
[47938]: https://github.com/pandas-dev/pandas/pull/47938
[51417]: https://github.com/pandas-dev/pandas/pull/51417
[28900]: https://github.com/pandas-dev/pandas/issues/28900
[35407]: https://github.com/pandas-dev/pandas/issues/35407
[53576]: https://github.com/pandas-dev/pandas/pull/53576
