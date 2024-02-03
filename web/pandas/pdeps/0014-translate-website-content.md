# PDEP-14: Publish translations of pandas.pydata.org

- Created: 01 February 2024
- Status: Under discussion
- Discussion: [#56301](https://github.com/pandas-dev/pandas/issues/56301)
              [#57204](https://github.com/pandas-dev/pandas/pull/57204)
- Author: [Albert Steppi](https://github.com/steppi),
- Revision: 1

## Abstract

The suggestion is to have official translations made for content of the core
project website [pandas.pydata.org](https://pandas.pydata.org) and offer
a low friction way for users to access these translations on the core
project website.

## Motivation, Scope, Usage, and Impact

There are many potential users with no or a low level of English proficiency
who could benefit from quality official translations of the Pandas website
content. Though translations for all documentation would be valuable,
producing and maintaining translations for such a large and oft-changing
collection of text would take an immense and sustained effort which may
be infeasible. The suggestion is instead to have translations made for only
a key set of pages from the core project website.

## Detailed Description and Implementation

The bulk of the work for setting up translation infrastructure, finding and
vetting translators, and working out how to publish translations, will fall
upon a cross-functional team funded by the [Scientific Python Community & Communications
Infrastructure grant](https://scientific-python.org/doc/scientific-python-community-and-communications-infrastructure-2022.pdf)
to work on adding translations for the main websites of all
[Scientific Python core projects](https://scientific-python.org/specs/core-projects/).
The hope is to minimize the burden on the core Pandas maintainers.

No translated content would be hosted within the Pandas repository itself.
Instead a separate GitHub repository could be set up containing the content
selected for translation. This repository could then be synced to the Crowdin
localization management platform as described in
[Crowdin's documentation](https://support.crowdin.com/github-integration/).
Crowdin would then provide a user interface for translators, and updates to
translations would be pushed to a feature branch, with completed translations
periodically merged into `main` after given approval by trusted
language-specific admin's working across the Scientific Python core projects
participating in the translation program. There will be no need for Pandas
maintainers to verify the quality of translations.

The result would be a repository containing parallel versions of content from
pandas.pydata.org, translated into various languages. Translated content could
then be pulled from this repository during generation of the Pandas website. A
low friction means of choosing between languages could then be added. Possibly a
drop-down language selector similar to what now exists for https://numpy.org, or
simple links similar to what now exists for https://www.sympy.org/en/index.html.
A developer supported by the "Scientific Python Community & Communications
Infrastructure grant" could assist with making the changes necessary for the
Pandas website to support publication of translations.

If desired, a cron job could be set up on the repository containing translated
content to check for relevant changes or updates to the Pandas website's content
and pull them if necessary. Translators could then receive a notification from
Crowdin that there are new strings to translate. This could help with the
process of keeping translations up to date.


### PDEP History

- 01 February 2024: Initial draft
