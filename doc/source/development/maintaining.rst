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

Triage is an important first step in addressing issues reported by the community, and even
partial contributions are a great way to help maintain pandas. Only remove the "Needs Triage"
tag once all of the steps below have been completed.

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
   new issue as "Duplicate" and close the issue with a link to the original issue.
   Make sure to still thank the reporter, and encourage them to chime in on the
   original issue, and perhaps try to fix it.

   If the new issue provides relevant information, such as a better or slightly
   different example, add it to the original issue as a comment or an edit to
   the original post.

4. **Is the issue minimal and reproducible**?

   For bug reports, we ask that the reporter provide a minimal reproducible
   example. See https://matthewrocklin.com/blog/work/2018/02/28/minimal-bug-reports
   for a good explanation. If the example is not reproducible, or if it's
   *clearly* not minimal, feel free to ask the reporter if they can provide
   an example or simplify the provided one. Do acknowledge that writing
   minimal reproducible examples is hard work. If the reporter is struggling,
   you can try to write one yourself and we'll edit the original post to include it.

   If a reproducible example can't be provided, add the "Needs info" label.

   If a reproducible example is provided, but you see a simplification,
   edit the original post with your simpler reproducible example.

   If this is a regression report, post the result of a ``git bisect`` run.
   More info on this can be found in the :ref:`maintaining.regressions` section.

   Ensure the issue exists on the main branch and that it has the "Needs Triage" tag
   until all steps have been completed. Add a comment to the issue once you have
   verified it exists on the main branch, so others know it has been confirmed.

5. **Is this a clearly defined feature request?**

   Generally, pandas prefers to discuss and design new features in issues, before
   a pull request is made. Encourage the submitter to include a proposed API
   for the new feature. Having them write a full docstring is a good way to
   pin down specifics.

   Tag new feature requests with "Needs Discussion", as we'll need a discussion
   from several pandas maintainers before deciding whether the proposal is in
   scope for pandas.

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

   If the issue is a regression report, add the "Regression" label and the next patch
   release milestone.

   Once you have completed the above, make sure to remove the "Needs Triage" label.

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
    git bisect run bash -c "python -m pip install -ve . --no-build-isolation -Ceditable-verbose=true; python t.py"

This finds the first commit that changed the behavior. The C extensions have to be
rebuilt at every step, so the search can take a while.

Exit bisect and rebuild the current version::

    git bisect reset
    python -m pip install -ve . --no-build-isolation -Ceditable-verbose=true

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
If you think an issue should be closed but are not completely sure, please apply
the "closing candidate" label and wait for other maintainers to take a look.

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
so passes with no response, thank them for their work and then either:

- close the pull request;
- push to the contributor's branch to push their work over the finish line (if
  you're part of ``pandas-core``). This can be helpful for pushing an important PR
  across the line, or for fixing a small merge conflict.

If closing the pull request, then please comment on the original issue that
"There's a stalled PR at #1234 that may be helpful.", and perhaps label the issue
as "Good first issue" if the PR was relatively close to being accepted.

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

If adding to ``pandas-core``, there are two additional steps:

3. Add the contributor to the pandas Google group.
4. Create a pull request to add the contributor's GitHub handle to ``pandas-dev/pandas/web/pandas/config.yml``.

The current list of core-team members is at
https://github.com/pandas-dev/pandas/blob/main/web/pandas/config.yml


.. _maintaining.merging:

Merging pull requests
---------------------

Only core team members can merge pull requests. We have a few guidelines.

1. You should typically not self-merge your own pull requests without approval.
   Exceptions include things like small changes to fix CI
   (e.g. pinning a package version). Self-merging with approval from other
   core team members is fine if the change is something you're very confident
   about.
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


.. _maintaining.release:

Release process
---------------

The release process makes a snapshot of pandas (a git commit) available to users with
a particular version number. After the release the new pandas version will be available
in the next places:

- Git repo with a `new tag <https://github.com/pandas-dev/pandas/tags>`_
- Source distribution in a `GitHub release <https://github.com/pandas-dev/pandas/releases>`_
- Pip packages in the `PyPI <https://pypi.org/project/pandas/>`_
- Conda packages in `conda-forge <https://anaconda.org/conda-forge/pandas>`_

The process for releasing a new version of pandas is detailed next section.

The instructions contain ``<version>`` which needs to be replaced with the version
to be released (e.g. ``1.5.2``). Also the branch to be released ``<branch>``, which
depends on whether the version being released is the release candidate of a new version,
or any other version. Release candidates are released from ``main``, while other
versions are released from their branch (e.g. ``1.5.x``).


Prerequisites
`````````````

In order to be able to release a new pandas version, the next permissions are needed:

- Merge rights to the `pandas <https://github.com/pandas-dev/pandas/>`_ and
  `pandas-feedstock <https://github.com/conda-forge/pandas-feedstock/>`_ repositories.
  For the latter, open a PR adding your GitHub username to the conda-forge recipe.
- Permissions to push to ``main`` in the pandas repository, to push the new tags.
- `Write permissions to PyPI <https://github.com/conda-forge/pandas-feedstock/pulls>`_.
- Access to our website / documentation server. Share your public key with the
  infrastructure committee to be added to the ``authorized_keys`` file of the main
  server user.
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

4. If not a release candidate, make sure all backporting pull requests to the
   branch being released are merged, and no merged pull requests are missing a
   backport (check the
   ["Still Needs Manual Backport"](https://github.com/pandas-dev/pandas/labels/Still%20Needs%20Manual%20Backport)
   label for this).

5. Create a new issue and milestone for the version after the one being released.
   If the release was a release candidate, we would usually want to create issues and
   milestones for both the next major/minor, and the next patch release. In the
   milestone of a patch release, we add the description ``on-merge: backport to <branch>``,
   so tagged PRs are automatically backported to the release branch by our bot.

6. Change the milestone of all issues and PRs in the milestone being released to the
   next milestone.

Release
```````

1. Create an empty commit and a tag in the last commit of the branch to be released::

    git checkout <branch>
    git pull --ff-only upstream <branch>
    git clean -xdf
    git commit --allow-empty --author="pandas Development Team <pandas-dev@python.org>" -m "RLS: <version>"
    git tag -a v<version> -m "Version <version>"  # NOTE that the tag is v1.5.2 with "v" not 1.5.2
    git push upstream <branch> --follow-tags

The docs for the new version will be built and published automatically with the docs job in the CI,
which will be triggered when the tag is pushed.

2. Only if the release is a release candidate, we want to create a new branch for it, immediately
   after creating the tag. For example, if we are releasing pandas 1.4.0rc0, we would like to
   create the branch 1.4.x to backport commits to the 1.4 versions. As well as create a tag to
   mark the start of the development of 1.5.0 (assuming it is the next version)::

    git checkout -b 1.4.x
    git push upstream 1.4.x
    git checkout main
    git commit --allow-empty -m "Start 1.5.0"
    git tag -a v1.5.0.dev0 -m "DEV: Start 1.5.0"
    git push upstream main --follow-tags

3. Download the source distribution and wheels from the `wheel staging area <https://anaconda.org/scientific-python-nightly-wheels/pandas>`_.
   Be careful to make sure that no wheels are missing (e.g. due to failed builds).

   Running scripts/download_wheels.sh with the version that you want to download wheels/the sdist for should do the trick.
   This script will make a ``dist`` folder inside your clone of pandas and put the downloaded wheels and sdist there::

    scripts/download_wheels.sh <VERSION>

   ATTENTION: this is currently not downloading *all* wheels, and you have to
   manually download the remainings wheels and sdist!

4. Create a `new GitHub release <https://github.com/pandas-dev/pandas/releases/new>`_:

   - Tag: ``<version>``
   - Title: ``pandas <version>``
   - Description: Copy the description of the last release of the same kind (release candidate, major/minor or patch release)
   - Files: ``pandas-<version>.tar.gz`` source distribution just generated
   - Set as a pre-release: Only check for a release candidate
   - Set as the latest release: Leave checked, unless releasing a patch release for an older version
     (e.g. releasing 1.4.5 after 1.5 has been released)

5. Upload wheels to PyPI::

    twine upload pandas/dist/pandas-<version>*.{whl,tar.gz} --skip-existing

6. The GitHub release will after some hours trigger an
   `automated conda-forge PR <https://github.com/conda-forge/pandas-feedstock/pulls>`_.
   (If you don't want to wait, you can open an issue titled ``@conda-forge-admin, please update version`` to trigger the bot.)
   Merge it once the CI is green, and it will generate the conda-forge packages.

   In case a manual PR needs to be done, the version, sha256 and build fields are the
   ones that usually need to be changed. If anything else in the recipe has changed since
   the last release, those changes should be available in ``ci/meta.yaml``.

Post-Release
````````````

1. Update symlinks to stable documentation by logging in to our web server, and
   editing ``/var/www/html/pandas-docs/stable`` to point to ``version/<X.Y>``
   for major and minor releases, or ``version/<X.Y>`` to ``version/<patch>`` for
   patch releases. The exact instructions are (replace the example version numbers by
   the appropriate ones for the version you are releasing):

   - Log in to the server and use the correct user.
   - ``cd /var/www/html/pandas-docs/``
   - For a major or minor release (assuming the ``/version/2.1.0/`` docs have been uploaded to the server):

     -  Create a new X.Y symlink to X.Y.Z: ``cd version; ln -sfn 2.1.0 2.1``
     -  Update stable symlink to point to X.Y: ``ln -sfn version/2.1 stable``

   - For a patch release (assuming the ``/version/2.1.3/`` docs have been uploaded to the server):

     - Update the X.Y symlink to the new X.Y.Z patch version: ``cd version; ln -sfn 2.1.3 2.1``
     - (the stable symlink should already be pointing to the correct X.Y version)

2. If releasing a major or minor release, open a PR in our source code to update
   ``web/pandas/versions.json``, to have the desired versions in the documentation
   dropdown menu.

3. Close the milestone and the issue for the released version.

4. Create a new issue for the next release, with the estimated date of release.

5. Open a PR with the placeholder for the release notes of the next version. See
   for example `the PR for 1.5.3 <https://github.com/pandas-dev/pandas/pull/49843/files>`_.
   Note that the template to use depends on whether it is a major, minor or patch release.

6. Announce the new release in the official channels (use previous announcements
   for reference):

   - The pandas-dev and pydata mailing lists
   - X, Mastodon, Telegram and LinkedIn

7. Update this release instructions to fix anything incorrect and to update about any
   change since the last release.

.. _governance documents: https://github.com/pandas-dev/pandas/blob/main/web/pandas/about/governance.md
.. _list of permissions: https://docs.github.com/en/organizations/managing-access-to-your-organizations-repositories/repository-roles-for-an-organization
