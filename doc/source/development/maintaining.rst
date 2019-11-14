.. _maintaining:

******************
Pandas Maintenance
******************

This guide is for pandas' maintainers. It may also be interesting to contributors
looking to understand the pandas development process and what steps are necessary
to become a maintainer.

Roles
-----

Pandas uses two levels of permissions: **triage** and **core** team members.

Triage members can label and close issues and pull requests.

Core team members can label and close issues and pull request, and can merge
pull requests.

Tasks
-----

Pandas is largely a volunteer project, so these tasks shouldn't be read as
"expectations" of triage and maintainers.

* Triage newly filed issues (see :ref:`maintaining.triage`)
* Review newly opened pull request
* Respond to updates on existing issues and pull requests

.. _maintaining.triage:

Issue Triage
------------

The issue tracker is many people's first interaction with the pandas project itself,
beyond just using the library. As such, we want it to be a welcoming, pleasant
experience.

Here's a typical workflow for triaging a newly opened issue.

1. **Is the necessary information provided?**

   Ideally reporters would fill out the issue template, but many don't.
   If crucial information (like the version of pandas they used), is missing
   feel free to ask for that and label the issue with "Needs info".

2. **Check for duplicates**

   We have many open issues. If a new issue is clearly a duplicate, label the
   new issue as "Duplicate" assign the milestone "No Action", and close the issue
   with a link to the original issue. Make sure to still thank the reporter, and
   encourage them to chime in on the original issue, and perhaps try to fix it.

   If the new issue provides relevant information, such as a better or slightly
   different example, add it to the original issue as a comment or an edit to
   the original post.

3. **Is the issue minimal and reproducible**?

   For bug reports, we ask that the reporter provide a minimal reproducible
   example. See http://matthewrocklin.com/blog/work/2018/02/28/minimal-bug-reports
   for a good explanation. If the example is not reproducible, or if it's
   *clearly* not minimal, feel free to ask the reporter if they can provide
   and example or simplify the provided one. Do acknowledge that writing
   minimal reproducible examples is hard work. If the reporter is struggling,
   you can try to write one yourself and we'll edit the original post to include it.

   If a reproducible example can't be provided, add the "Needs info" label.

4. **Feature Requests**

   Generally, pandas prefers to discuss and design new features in issues, before
   a pull request is made. Encourage the submitter to include a proposed API
   for the new feature. Having them write a full docstring is a good exercise.

5. **Usage Questions**

   We prefer that usage questions are asked on StackOverflow. If it's easy to
   answer, feel free to link to the relevant documentation section, let them
   know that in the future this kind of question should be on StackOverflow,
   and close the issue.

6. **Labels and Milestones**

   Apply the relevant labels. This is a bit of an art, and comes with experience.
   Look at similar issues to get a feel for how things are labeled.

   If there issue is clearly defined and the fix seems relatively straightforward,
   label the issue as "Good first issue".

   Typically, new issues will be assigned the "Contributions welcome" milestone,
   unless it's know that this issue should be addressed in a specific release (say
   because it's a large regression).

Closing Issues
--------------

Be delicate here: many people interpret closing an issue as us saying that the
conversation is over. It's typically best to give the reporter some time to
self-close their issue if it's determined that the behavior is not a bug,
or the feature is out of scope. Sometimes reporters just go away though, and
we'll close the issue after the conversation has died.

Reviewing Pull Requests
-----------------------

Anybody can review a pull request: regular contributors, triagers, or core-team
members. Here are some guidelines to check.

* Tests should be in a sensible location.
* New public APIs should be included in ``doc/source/reference/``.
* New / changed API should use the ``versionadded`` or ``versionchanged`` directives in the docstring.
* User-facing changes should have a whatsnew in the appropriate file.
* Regression tests should reference the original GitHub issue number like ``# GH-1234``.

.. _people: https://github.com/pandas-dev/pandas-governance/blob/master/people.md
