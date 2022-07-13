# Project governance

The official version of this document, along with a list of
individuals and institutions in the roles defined in the governance
section below, is contained in the
[Project governance](https://pandas.pydata.org/about/governance.html)
page of the pandas website.

## The Project

The pandas Project (The Project) is an open source software project affiliated
with the 501(c)3 NumFOCUS Foundation. The goal of The Project is to develop open
source software for data ingest, data preparation, data analysis, and data
visualization for the Python programming language. The Software developed by
The Project is released under the BSD (or similar) open source license,
developed openly and hosted in public GitHub repositories under the [pandas
GitHub organization](https://github.com/pandas-dev). Examples of Project Software
include the main pandas code repository and the pandas-stubs library.

Through its affiliation with NumFOCUS, The Project has the right to receive
tax-deductible donations in the United States of America.

The Project is developed by a team of distributed developers, called
Contributors. Contributors are individuals who have contributed code,
documentation, designs or other work to one or more Project repositories.
Anyone can be a Contributor. Contributors can be affiliated with any legal
entity or none. Contributors participate in the project by submitting,
reviewing and discussing GitHub Pull Requests and Issues and participating in
open and public Project discussions on GitHub, mailing lists, and
elsewhere. The foundation of Project participation is openness and
transparency.

Here is a list of the current Contributors to the main pandas repository:

[https://github.com/pandas-dev/pandas/graphs/contributors](https://github.com/pandas-dev/pandas/graphs/contributors)

There are also many other Contributors listed in the logs of other repositories of
the pandas project.

The Project Community consists of all Contributors and Users of the Project.
Contributors work on behalf of and are responsible to the larger Project
Community and we strive to keep the barrier between Contributors and Users as
low as possible.

The Project is formally affiliated with the 501(c)3 NumFOCUS Foundation
([https://numfocus.org](https://numfocus.org)), which serves as its fiscal
sponsor, may hold project trademarks and other intellectual property, helps
manage project donations and acts as a parent legal entity. NumFOCUS is the
only legal entity that has a formal relationship with the project (see
Institutional Partners section below).

## Community members

The pandas community is composed of a diverse group of stakeholders, such as
developers, contributors, individual and corporate users, etc. There are
some groups which have specific responsibilities. We list them next.

### Active maintainers

Active maintainers (aka the core developer team) are contributors of the
project who made significant contributions in the form of code, reviews,
software design, documentation etc. Their role in the project is to
advance the pandas software and goals.

**Membership**: Contributors to the pandas project become maintainers after
showing significant contributions over a period of one year, and with
consensus from the rest of activate maintainers

See the list of active maintainers [here](team.html#maintainers).

### BDFL (Benevolent Dictator For Life)

The figure of the BDFL exist to be able to unblock situations where a decision
needs to be made, and consensus or voting has failed. In such situations, the
BDFL will make the final decision.

**Membership**: Wes McKinney, as original creator of pandas has been the BDFL
of the project. In the event of Wes stepping down as BDFL, maintainers will
make a decision about whether to appoint a new BDFL or change the governance
type.

### Finances committee

The role of the members of the finances committee is to approve the spending
of pandas funds. Decisions in general will be made together with the rest of
active maintainers. Committee members will be responsible to make the final
decisions, and formally approve payments.

**Membership**: The committee will have 5 members, who will be selected by active
maintainers. Some constraints exists regarding committee membership:

- Members must be active maintainers
- No more than two committee members can be employed directly or indirectly
  by the same employer
- Committee members should not have conflicts of interest that could prevent
  them to make the best decisions in the interest of the project. This includes
  maintainers who receive significant payments from pandas funds

### Code of conduct committee

The role of the committee is to make sure pandas is as open, transparent and
inclusive as it aims to be by its values. In particular, the committee will
monitor and respond to any possible violation of our code of conduct. And
will publish regular summaries about violation reports.

**Membership**: Any members of the community can be part of the committee.
The committee will have 5 members, who will be selected by active maintainers.
The next constraints must be satisfied:

- The committee should aim to be as diverse as reasonably possible, to be able
  to make decisions based on a variety of points of views. In particular, the
  committee should not have more than 3 members of the same gender, or more
  than two members from the same geography (continent). Ideally the committee
  will also be diverse in other ways such as religion, political views,
  age, etc.
- No more than two members of the committee should be pandas maintainers.

### Inactive maintainers

Inactive maintainers are former active maintainers who are not interested or
not available to continue contributing to pandas in a regular way. If they
decide to participate in a discussion, they will still be considered active
maintainers for that discussion, but otherwise they are not expected to be part
of the decision making of the project, not have commit rights to the pandas
repositories, or be in the maintainers distribution list.

**Membership**: Active maintainers become inactive by their own decision.
Inactive maintainers can become active again if they are interested.

### NumFOCUS

[NumFOCUS](https://numfocus.org) is the fiscal sponsor of the pandas project.
As such, NumFOCUS is the legal and financial entity of the project, being the
owner of pandas trademarks and copyrights, and the legal entity of the
project for financial and tax reasons. NumFOCUS also helps promote pandas, and
find synergies with other projects of the ecosystem.

### Sponsors

Sponsors are institutions (companies, non-profits, universities, government
agencies, etc) that contribute to the pandas project. The main types of
sponsors are institutions employing people who work in pandas as part of their
job. And institutions funding the project. Sponsors will have advantages like
being listed in the pandas website, being mentioned in pandas channels
such as the blog or social media, or having direct communication with the
pandas maintainers, other than the usual channels. And others agreed by
active maintainers.

**Membership**: Institutions become sponsors if they employ a person to work
on pandas at least one day per week. Or if they provide funds to the project
(in money or in kind) of value of at least $10,000. Institutions stop
being considered sponsors after one year since the last action that made them
sponsors.

## Private communications of the Core Team

Unless specifically required, all Core Team discussions and activities will be
public and done in collaboration and discussion with the Project Contributors
and Community. The Core Team will have a private mailing list that will be used
sparingly and only when a specific matter requires privacy. When private
communications and decisions are needed, the Core Team will do its best to
summarize those to the Community after eliding personal/private/sensitive
information that should not be posted to the public internet.

## Breach

Non-compliance with the terms of the governance documents shall be reported to
the Core Team either through public or private channels as deemed appropriate.

## Changing the Governance

Changes to the governance are submitted via a GitHub pull request to The Project's
[governance page](https://github.com/pandas-dev/pandas/blob/main/web/pandas/about/governance.md).
The pull request is then refined in response to public comment and review, with
the goal being consensus in the community.  After this open period, a Core Team
Member proposes to the Core Team that the changes be ratified and the pull
request merged (accepting the proposed changes) or proposes that the pull
request be closed without merging (rejecting the proposed changes). The Member
should state the final commit hash in the pull request being proposed for
acceptance or rejection and briefly summarize the pull request. A minimum of
80% of the Core Team must vote and at least 2/3 of the votes must be positive
to carry out the proposed action (fractions of a vote rounded up to the nearest
integer). Since the BDFL holds ultimate authority in The Project, the BDFL has
authority to act alone in accepting or rejecting changes or overriding Core
Team decisions.
