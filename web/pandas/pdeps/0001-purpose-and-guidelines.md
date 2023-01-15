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

Bug fixes and conceptually minor changes (e.g. adding a parameter to a function)
are out of the scope of PDEPs. A PDEP should be used for changes that are not
immediate and not obvious, and are expected to require a significant amount of
discussion and require detailed documentation before being implemented.

PDEP are appropriate for user facing changes, internal changes and organizational
discussions. Examples of topics worth a PDEP could include moving a module from
pandas to a separate repository, a refactoring of the pandas block manager or
a proposal of a new code of conduct.

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

Anyone can propose a PDEP, but in most cases developers of pandas itself and related
projects are expected to author PDEPs. If you are unsure if you should be opening
an issue or creating a PDEP, it's probably safe to start by
[opening an issue](https://github.com/pandas-dev/pandas/issues/new/choose), which can
be eventually moved to a PDEP.

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

The initial status of a PDEP will be `Status: Under discussion`. This will be changed
to `Status: Accepted` when the PDEP is ready and have the approval of the core team.

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

### PDEP-1 History

- 3 August 2022: Initial version
