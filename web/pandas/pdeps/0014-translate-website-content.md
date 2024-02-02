# PDEP-14: Publish translations of pandas.pydata.org

- Created: 01 February 2024
- Status: Under discussion
- Discussion: [#56301](https://github.com/pandas-dev/pandas/issues/56301)
              [#57204](https://github.com/pandas-dev/pandas/pull/57204)
- Author: [Albert Steppi](https://github.com/steppi),
- Revision: 1

## Abstract

The suggestion is to have official translations made for content of the core
project website [pandas.pydata.org](https://pandas.pydata.org) and provide a
language drop-down selector on [pandas.pydata.org](https://pandas.pydata.org)
similar to what currently exists at [numpy.org](https://numpy.org).


## Motivation and Scope

Pandas is a foundational package in the Scientific Python ecosystem and there
are many potential users with no or low English proficiency who would benefit
from having high quality information about Pandas available in their native
language.

Translation of all content presents considerable challenge due to its sheer
volume and due to the tendency for technical documentation to exist in a state
of flux. The suggestion is to have translations for a targeted subset, selected:

- from things which are relatively stable to reduce the ongoing burden of
  keeping translations up to date.
- to maximize the benefit to users and potential users who currently have no or
  a low level of English proficiency, given the person-hours and resources that
  are likely to be available now and into the future.

Consideration of what subset of content would be most useful for users with
no or a low level of English proficiency could be a guiding principal to help
select what information should be available on the core project website, outside
of the technical documentation.

## Detailed Description

The following is a list of all pages on the core project website which are sourced
from markdown files at https://github.com/pandas-dev/pandas/tree/main/web/pandas.

- Landing page: https://pandas.pydata.org
- About pandas: https://pandas.pydata.org/about
- Project roadmap: https://pandas.pydata.org/about/roadmap.html
- Governance: https://pandas.pydata.org/about/governance.html
- Team: https://pandas.pydata.org/about/team.html
- Sponsors: https://pandas.pydata.org/about/sponsors.html
- Citing and logo: https://pandas.pydata.org/about/citing.html
- Getting started: https://pandas.pydata.org/getting_started.html
- Code of conduct: https://pandas.pydata.org/community/coc.html
- Ecosystem: https://pandas.pydata.org/community/ecosystem.html
- Contribute: https://pandas.pydata.org/contribute.html

Provisionally, the suggestion is for all of this content to be translated with
the possible exception of the "Project roadmap", which may be of limited
interest to new users.  Currently the "Getting started" section may be of
limited utility to users unable to engage with the externally linked content. In
the "Project roadmap" within the subsection labeled "Documentation improvements"
there is a stated goal to:

*Improve the "Getting Started" documentation, designing and writing learning
 paths for users different backgrounds (e.g. brand new to programming, familiar
 with other languages like R, already familiar with Python).*

It is recommended that this goal be accomplished alongside translation work in
order to make this page more useful to those with no or low English proficiency.
This would also prevent the need for retranslation if this goal were to be
accomplished after the original translation work is completed.

A language selection drop-down should be added to the navigation-bar similar to
what exists at https://numpy.org.


## Usage and Impact

The primary impact would be lowering the barrier to entry for non-English
speakers to get started using Pandas and moving along the path towards learning
to use it skillfully.

In 2022 it was estimated that there were approximately 400 million native
speakers of English and between 1.5 - 2 billion people who speak English as a
second language worldwide
[Wikipedia](https://web.archive.org/web/20240129080609/https://en.wikipedia.org/wiki/English-speaking_world).
With an estimated world population of over 8 billion people, this leaves many
for whom the Pandas core website is not directly accessible. Pandas is an
important piece of software infrastructure for data manipulation and analysis
with utility beyond the English speaking world. There is a vast population of
users and potential users who could benefit from having official information
about Pandas published in their native language.

Although automated translation tools can help those with no or low English
proficiency access the content of the Pandas website, these tools often still
struggle with the technical and jargon-laden language of scientific
software. This was evinced during the translation of https://numpy.org.
Automatic translation tools are invaluable as a starting point for human
translators, but human translators remain important to ensure accuracy.

## Implementation

The bulk of the work for setting up translation infrastructure, finding and
vetting translators, and working out how to publish translations, will fall
upon a cross-functional team funded by the [Scientific Python Community & Communications
Infrastructure grant](https://scientific-python.org/doc/scientific-python-community-and-communications-infrastructure-2022.pdf)
to work on adding translations for the main websites of all
[Scientific Python core projects](https://scientific-python.org/specs/core-projects/).
The goal is to minimize the burden on the core Pandas maintainers.

A GitHub repository should be set up to mirror content from the core webpage
which is selected for translation. A GitHub action should be set up to keep
the mirrored repository up-to-date. Either an action within the main Pandas
repo which pushes updates to the mirror, or a cron in the mirror which polls
for relevant updates in Pandas repo and pulls them when necessary.

The mirrored repository would then be synced to the Crowdin localization
management platform as described in
[Crowdin's documentation](https://support.crowdin.com/github-integration/).
There would be separate folders within the mirror repository, one for each target
language, with the content initially untranslated.
Crowdin would then provide a user interface for translators, and updates
to translations would be pushed to the branch `l10n_main` on the mirrored
repository. Periodically, manual pull requests would be made to the main Pandas
repo, adding translated content within folders alongside of the English content.

Translations will be managed within an enterprise Crowdin organization created for
Scientific Python localization projects. Access to this organization is
invite-only, and translators will be vetted to help safe-guard against the
spamming of low quality or inflammatory translations. Approval from a trusted
admin would be required before translations are merged into the main Pandas
repo.

A language drop-down selector will need to be added to the navigation-bar of
the Pandas website. The plan is for development of a generic solution that
can be reused for all Scientific Python website translations.


### PDEP History

- 01 February 2024: Initial draft
