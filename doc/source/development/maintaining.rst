.. _maintaining:

******************
Pandas Maintenance
******************

This guide is for pandas' maintainers. It may also be interesting to contributors
looking to understand the pandas development process and what steps are necessary
to become a maintainer.

The main contributing guide is available at :ref:`contributing`.

Roles
-----

Pandas uses two levels of permissions: **triage** and **core** team members.

Triage members can label and close issues and pull requests.

Core team members can label and close issues and pull request, and can merge
pull requests.

GitHub publishes the full `list of permissions`_.

Tasks
-----

Pandas is largely a volunteer project, so these tasks shouldn't be read as
"expectations" of triage and maintainers. Rather, they're general descriptions
of what it means to be a maintainer.

* Triage newly filed issues (see :ref:`maintaining.triage`)
* Review newly opened pull requests
* Respond to updates on existing issues and pull requests
* Drive discussion and decisions on stalled issues and pull requests
* Provide experience / wisdom on API design questions to ensure consistency and maintainability
* Project organization (run / attend developer meetings, represent pandas)

https://matthewrocklin.com/blog/2019/05/18/maintainer may be interesting background
reading.

.. _maintaining.triage:

Issue Triage
------------


Here's a typical workflow for triaging a newly opened issue.

1. **Thank the reporter for opening an issue**

   The issue tracker is many people's first interaction with the pandas project itself,
   beyond just using the library. As such, we want it to be a welcoming, pleasant
   experience.

2. **Is the necessary information provided?**

   Ideally reporters would fill out the issue template, but many don't.
   If crucial information (like the version of pandas they used), is missing
   feel free to ask for that and label the issue with "Needs info". The
   report should follow the guidelines in :ref:`contributing.bug_reports`.
   You may want to link to that if they didn't follow the template.

   Make sure that the title accurately reflects the issue. Edit it yourself
   if it's not clear.

3. **Is this a duplicate issue?**

   We have many open issues. If a new issue is clearly a duplicate, label the
   new issue as "Duplicate" assign the milestone "No Action", and close the issue
   with a link to the original issue. Make sure to still thank the reporter, and
   encourage them to chime in on the original issue, and perhaps try to fix it.

   If the new issue provides relevant information, such as a better or slightly
   different example, add it to the original issue as a comment or an edit to
   the original post.

4. **Is the issue minimal and reproducible**?

   For bug reports, we ask that the reporter provide a minimal reproducible
   example. See https://matthewrocklin.com/blog/work/2018/02/28/minimal-bug-reports
   for a good explanation. If the example is not reproducible, or if it's
   *clearly* not minimal, feel free to ask the reporter if they can provide
   and example or simplify the provided one. Do acknowledge that writing
   minimal reproducible examples is hard work. If the reporter is struggling,
   you can try to write one yourself and we'll edit the original post to include it.

   If a reproducible example can't be provided, add the "Needs info" label.

   If a reproducible example is provided, but you see a simplification,
   edit the original post with your simpler reproducible example.

5. **Is this a clearly defined feature request?**

   Generally, pandas prefers to discuss and design new features in issues, before
   a pull request is made. Encourage the submitter to include a proposed API
   for the new feature. Having them write a full docstring is a good way to
   pin down specifics.

   We'll need a discussion from several pandas maintainers before deciding whether
   the proposal is in scope for pandas.

6. **Is this a usage question?**

   We prefer that usage questions are asked on StackOverflow with the pandas
   tag. https://stackoverflow.com/questions/tagged/pandas

   If it's easy to answer, feel free to link to the relevant documentation section,
   let them know that in the future this kind of question should be on
   StackOverflow, and close the issue.

7. **What labels and milestones should I add?**

   Apply the relevant labels. This is a bit of an art, and comes with experience.
   Look at similar issues to get a feel for how things are labeled.

   If the issue is clearly defined and the fix seems relatively straightforward,
   label the issue as "Good first issue".

   Typically, new issues will be assigned the "Contributions welcome" milestone,
   unless it's know that this issue should be addressed in a specific release (say
   because it's a large regression).

.. _maintaining.closing:

Closing Issues
--------------

Be delicate here: many people interpret closing an issue as us saying that the
conversation is over. It's typically best to give the reporter some time to
respond or self-close their issue if it's determined that the behavior is not a bug,
or the feature is out of scope. Sometimes reporters just go away though, and
we'll close the issue after the conversation has died.

Reviewing Pull Requests
-----------------------

Anybody can review a pull request: regular contributors, triagers, or core-team
members. Here are some guidelines to check.

* Tests should be in a sensible location.
* New public APIs should be included somewhere in ``doc/source/reference/``.
* New / changed API should use the ``versionadded`` or ``versionchanged`` directives in the docstring.
* User-facing changes should have a whatsnew in the appropriate file.
* Regression tests should reference the original GitHub issue number like ``# GH-1234``.

Cleaning up old Issues
----------------------

Every open issue in pandas has a cost. Open issues make finding duplicates harder,
and can make it harder to know what needs to be done in pandas. That said, closing
issues isn't a goal on its own. Our goal is to make pandas the best it can be,
and that's best done by ensuring that the quality of our open issues is high.

Occasionally, bugs are fixed but the issue isn't linked to in the Pull Request.
In these cases, comment that "This has been fixed, but could use a test." and
label the issue as "Good First Issue" and "Needs Test".

If an older issue doesn't follow our issue template, edit the original post to
include a minimal example, the actual output, and the expected output. Uniformity
in issue reports is valuable.

If an older issue lacks a reproducible example, label it as "Needs Info" and
ask them to provide one (or write one yourself if possible). If one isn't
provide reasonably soon, close it according to the policies in :ref:`maintaining.closing`.

Cleaning up old Pull Requests
-----------------------------

Occasionally, contributors are unable to finish off a pull request.
If some time has passed (two weeks, say) since the last review requesting changes,
gently ask if they're still interested in working on this. If another two weeks or
so passes with no response, thank them for their work and close the pull request.
Comment on the original issue that "There's a stalled PR at #1234 that may be
helpful.", and perhaps label the issue as "Good first issue" if the PR was relatively
close to being accepted.

Additionally, core-team members can push to contributors branches. This can be
helpful for pushing an important PR across the line, or for fixing a small
merge conflict.

Becoming a pandas maintainer
----------------------------

The full process is outlined in our `governance documents`_. In summary,
we're happy to give triage permissions to anyone who shows interest by
being helpful on the issue tracker.

The current list of core-team members is at
https://github.com/pandas-dev/pandas-governance/blob/master/people.md

.. _governance documents: https://github.com/pandas-dev/pandas-governance
.. _list of permissions: https://help.github.com/en/github/setting-up-and-managing-organizations-and-teams/repository-permission-levels-for-an-organization