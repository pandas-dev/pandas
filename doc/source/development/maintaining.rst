.. _maintaining:

******************
pandas maintenance
******************

This guide is for pandas' maintainers. It may also be interesting to contributors
looking to understand the pandas development process and what steps are necessary
to become a maintainer.

The main contributing guide is available at :ref:`contributing`.

Roles
-----

pandas uses two levels of permissions: **triage** and **core** team members.

Triage members can label and close issues and pull requests.

Core team members can label and close issues and pull request, and can merge
pull requests.

GitHub publishes the full `list of permissions`_.

Tasks
-----

pandas is largely a volunteer project, so these tasks shouldn't be read as
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

Issue triage
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

.. _maintaining.regressions:

Investigating regressions
-------------------------

Regressions are bugs that unintentionally break previously working code. The common way
to  investigate regressions is by using
`git bisect <https://git-scm.com/docs/git-bisect>`_,
which finds the first commit that introduced the bug.

For example: a user reports that ``pd.Series([1, 1]).sum()`` returns ``3``
in pandas version ``1.5.0`` while in version ``1.4.0`` it returned ``2``. To begin,
create a file ``t.py`` in your pandas directory, which contains

.. code-block:: python

    import pandas as pd
    assert pd.Series([1, 1]).sum() == 2

and then run::

    git bisect start
    git bisect good v1.4.0
    git bisect bad v1.5.0
    git bisect run bash -c "python setup.py build_ext -j 4; python t.py"

This finds the first commit that changed the behavior. The C extensions have to be
rebuilt at every step, so the search can take a while.

Exit bisect and rebuild the current version::

    git bisect reset
    python setup.py build_ext -j 4

Report your findings under the corresponding issue and ping the commit author to get
their input.

.. note::
    In the ``bisect run`` command above, commits are considered good if ``t.py`` exits
    with ``0`` and bad otherwise. When raising an exception is the desired behavior,
    wrap the code in an appropriate ``try/except`` statement. See :issue:`35685` for
    more examples.

.. _maintaining.closing:

Closing issues
--------------

Be delicate here: many people interpret closing an issue as us saying that the
conversation is over. It's typically best to give the reporter some time to
respond or self-close their issue if it's determined that the behavior is not a bug,
or the feature is out of scope. Sometimes reporters just go away though, and
we'll close the issue after the conversation has died.

.. _maintaining.reviewing:

Reviewing pull requests
-----------------------

Anybody can review a pull request: regular contributors, triagers, or core-team
members. But only core-team members can merge pull requests when they're ready.

Here are some things to check when reviewing a pull request.

* Tests should be in a sensible location: in the same file as closely related tests.
* New public APIs should be included somewhere in ``doc/source/reference/``.
* New / changed API should use the ``versionadded`` or ``versionchanged`` directives in the docstring.
* User-facing changes should have a whatsnew in the appropriate file.
* Regression tests should reference the original GitHub issue number like ``# GH-1234``.
* The pull request should be labeled and assigned the appropriate milestone (the next patch release
  for regression fixes and small bug fixes, the next minor milestone otherwise)
* Changes should comply with our :ref:`policies.version`.


.. _maintaining.backporting:

Backporting
-----------

pandas supports point releases (e.g. ``1.4.3``) that aim to:

1. Fix bugs in new features introduced in the first minor version release.

  * e.g. If a new feature was added in ``1.4`` and contains a bug, a fix can be applied in ``1.4.3``

2. Fix bugs that used to work in a few minor releases prior. There should be agreement between core team members that a backport is appropriate.

  * e.g. If a feature worked in ``1.2`` and stopped working since ``1.3``, a fix can be applied in ``1.4.3``.

Since pandas minor releases are based on GitHub branches (e.g. point release of ``1.4`` are based off the ``1.4.x`` branch),
"backporting" means merging a pull request fix to the ``main`` branch and correct minor branch associated with the next point release.

By default, if a pull request is assigned to the next point release milestone within the GitHub interface,
the backporting process should happen automatically by the ``@meeseeksdev`` bot once the pull request is merged.
A new pull request will be made backporting the pull request to the correct version branch.
Sometimes due to merge conflicts, a manual pull request will need to be made addressing the code conflict.

If the bot does not automatically start the backporting process, you can also write a GitHub comment in the merged pull request
to trigger the backport::

    @meeseeksdev backport version-branch

This will trigger a workflow which will backport a given change to a branch
(e.g. @meeseeksdev backport 1.4.x)

Cleaning up old issues
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

Cleaning up old pull requests
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

The required steps for adding a maintainer are:

1. Contact the contributor and ask their interest to join.
2. Add the contributor to the appropriate `GitHub Team <https://github.com/orgs/pandas-dev/teams>`_ if accepted the invitation.

  * ``pandas-core`` is for core team members
  * ``pandas-triage`` is for pandas triage members

3. Add the contributor to the pandas Google group.
4. Create a pull request to add the contributor's GitHub handle to ``pandas-dev/pandas/web/pandas/config.yml``.
5. Create a pull request to add the contributor's name/GitHub handle to the `governance document <https://github.com/pandas-dev/pandas-governance/blob/master/people.md>`_.

The current list of core-team members is at
https://github.com/pandas-dev/pandas-governance/blob/master/people.md


.. _maintaining.merging:

Merging pull requests
---------------------

Only core team members can merge pull requests. We have a few guidelines.

1. You should typically not self-merge your own pull requests. Exceptions include
   things like small changes to fix CI (e.g. pinning a package version).
2. You should not merge pull requests that have an active discussion, or pull
   requests that has any ``-1`` votes from a core maintainer. pandas operates
   by consensus.
3. For larger changes, it's good to have a +1 from at least two core team members.

In addition to the items listed in :ref:`maintaining.closing`, you should verify
that the pull request is assigned the correct milestone.

Pull requests merged with a patch-release milestone will typically be backported
by our bot. Verify that the bot noticed the merge (it will leave a comment within
a minute typically). If a manual backport is needed please do that, and remove
the "Needs backport" label once you've done it manually. If you forget to assign
a milestone before tagging, you can request the bot to backport it with:

.. code-block:: console

   @Meeseeksdev backport <branch>


.. _maintaining.asv-machine:

Benchmark machine
-----------------

The team currently owns dedicated hardware for hosting a website for pandas' ASV performance benchmark. The results
are published to https://asv-runner.github.io/asv-collection/pandas/

Configuration
`````````````

The machine can be configured with the `Ansible <http://docs.ansible.com/ansible/latest/index.html>`_ playbook in https://github.com/tomaugspurger/asv-runner.

Publishing
``````````

The results are published to another GitHub repository, https://github.com/tomaugspurger/asv-collection.
Finally, we have a cron job on our docs server to pull from https://github.com/tomaugspurger/asv-collection, to serve them from ``/speed``.
Ask Tom or Joris for access to the webserver.

Debugging
`````````

The benchmarks are scheduled by Airflow. It has a dashboard for viewing and debugging the results. You'll need to setup an SSH tunnel to view them

    ssh -L 8080:localhost:8080 pandas@panda.likescandy.com


.. _maintaining.release:

Release process
---------------

The release process makes a snapshot of pandas (a git commit) available to users with
a particular version number. After the release the new pandas version will be available
in the next places:

- Git repo with a [new tag](https://github.com/pandas-dev/pandas/tags)
- Source distribution in a [GitHub release](https://github.com/pandas-dev/pandas/releases)
- Pip packages in the [PyPI](https://pypi.org/project/pandas/)
- Conda/Mamba packages in [conda-forge](https://anaconda.org/conda-forge/pandas)

The process for releasing a new version of pandas is detailed next section.

The instructions contain ``<version>`` which needs to be replaced with the version
to be released (e.g. ``1.5.2``). Also the branch to be released ``<branch>``, which
depends on whether the version being released is the release candidate of a new version,
or any other version. Release candidates are released from ``main``, while other
versions are released from their branch (e.g. ``1.5.x``).


Prerequisites
`````````````

In order to be able to release a new pandas version, the next permissions are needed:

- Merge rights to the [pandas](https://github.com/pandas-dev/pandas/),
  [pandas-wheels](https://github.com/MacPython/pandas-wheels), and
  [pandas-feedstock](https://github.com/conda-forge/pandas-feedstock/) repositories.
- Permissions to push to main in the pandas repository, to push the new tags.
- Write permissions to [PyPI](https://github.com/conda-forge/pandas-feedstock/pulls)
- Access to the social media accounts, to publish the announcements.

Pre-release
```````````

1. Agree with the core team on the next topics:

   - Release date (major/minor releases happen usually every 6 months, and patch releases
     monthly until x.x.5, just before the next major/minor)
   - Blockers (issues and PRs that must be part of the release)
   - Next version after the one being released

2. Update and clean release notes for the version to be released, including:

   - Set the final date of the release
   - Remove any unused bullet point
   - Make sure there are no formatting issues, typos, etc.

3. Make sure the CI is green for the last commit of the branch being released.

4. If not a release candidate, make sure all backporting pull requests to the branch
   being released are merged.

5. Create a new issue and milestone for the version after the one being released.
   If the release was a release candidate, we would usually want to create issues and
   milestones for both the next major/minor, and the next patch release. In the
   milestone of a patch release, we add the description ``on-merge: backport to <branch>``,
   so tagged PRs are automatically backported to the release branch by our bot.

6. Change the milestone of all issues and PRs in the milestone being released to the
   next milestone.

Release
```````

1. Create an empty commit and a tag in the last commit of the branch to be released:

    git checkout <branch>
    git pull --ff-only upstream <branch>
    git clean -xdf
    git commit --allow-empty --author="Pandas Development Team <pandas-dev@python.org>" -m "RLS: <version>"
    git tag -a v<version> -m "Version <version>"  # NOTE that the tag is v1.5.2 with "v" not 1.5.2
    git push upstream <branch> --follow-tags

The docs for the new version will be built and published automatically with the docs job in the CI,
which will be triggered when the tag is pushed.

2. Only if the release is a release candidate, we want to create a new branch for it, immediately
   after creating the tag. For example, if we are releasing pandas 1.4.0rc0, we would like to
   create the branch 1.4.x to backport commits to the 1.4 versions. As well as create a tag to
   mark the start of the development of 1.5.0 (assuming it is the next version):

    git checkout -b 1.4.x
    git push upstream 1.4.x
    git checkout main
    git commit --allow-empty -m "Start 1.5.0"
    git tag -a v1.5.0.dev0 -m "DEV: Start 1.5.0"
    git push upstream main --follow-tags

3. Build the source distribution (git must be in the tag commit):

    ./setup.py sdist --formats=gztar --quiet

4. Create a [new GitHub release](https://github.com/pandas-dev/pandas/releases/new):

   - Title: ``Pandas <version>``
   - Tag: ``<version>``
   - Files: ``pandas-<version>.tar.gz`` source distribution just generated
   - Description: Copy the description of the last release of the same kind (release candidate, major/minor or patch release)
   - Set as a pre-release: Only check for a release candidate
   - Set as the latest release: Leave checked, unless releasing a patch release for an older version
     (e.g. releasing 1.4.5 after 1.5 has been released)

5. The GitHub release will after some hours trigger an
   [automated conda-forge PR](https://github.com/conda-forge/pandas-feedstock/pulls).
   Merge it once the CI is green, and it will generate the conda-forge packages.

6. Packages for supported versions in PyPI are built in the
   [MacPython repo](https://github.com/MacPython/pandas-wheels).
   Open a PR updating the build commit to the released version, and merge it once the
   CI is green.

    git checkout master
    git pull --ff-only upstream master
    git checkout -B RLS-<version>
    sed -i 's/BUILD_COMMIT: "v.*/BUILD_COMMIT: "'v<version>'"/' azure/windows.yml azure/posix.yml
    sed -i 's/BUILD_COMMIT="v.*/BUILD_COMMIT="'v<version>'"/' .travis.yml
    git commit -am "RLS <version>"
    git push -u origin RLS-<version>

7. Download all wheels from the Anaconda repository where MacPython uploads them:
   https://anaconda.org/multibuild-wheels-staging/pandas/files?version=<version>
   to the ``dist/`` directory in the local pandas copy.

8. Upload wheels to PyPI:

    twine upload pandas/dist/pandas-<version>*.{whl,tar.gz} --skip-existing

Post-Release
````````````

1. Update symlink to stable documentation by logging in to our web server, and
   editing ``/var/www/html/pandas-docs/stable`` to point to ``version/<latest-version>``.

2. If releasing a major or minor release, open a PR in our source code to update
   ``web/pandas/versions.json``, to have the desired versions in the documentation
   dropdown menu.

3. Close the milestone and the issue for the released version.

4. Create a new issue for the next release, with the estimated date of release.

5. Open a PR with the placeholder for the release notes of the next version. See
   for example [the PR for 1.5.3](https://github.com/pandas-dev/pandas/pull/49843/files).

6. Announce the new release in the official channels (use previous announcements
   for reference):

    - The pandas-dev and pydata mailing lists
    - Twitter, Mastodon and Telegram


.. _governance documents: https://github.com/pandas-dev/pandas-governance
.. _list of permissions: https://docs.github.com/en/organizations/managing-access-to-your-organizations-repositories/repository-roles-for-an-organization
