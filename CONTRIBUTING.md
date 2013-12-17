###Guidelines

All contributions, bug reports, bug fixes, documentation improvements,
enhancements and ideas are welcome.

The [GitHub "issues" tab](https://github.com/pydata/pandas/issues)
contains some issues labeled "Good as first PR"; Look those up if you're
looking for a quick way to help out.

#### Bug Reports

  - Please include a short, self-contained Python snippet reproducing the problem.
  You can have the code formatted nicely by using [GitHub Flavored Markdown](http://github.github.com/github-flavored-markdown/) :

        ```python

        print("I ♥ pandas!")

        ```

  - Specify the pandas version used and those of it's dependencies. You can simply include   the output of
    [`ci/print_versions.py`](https://github.com/pydata/pandas/blob/master/ci/print_versions.py).
  - Explain what the expected behavior was, and what you saw instead.

#### Pull Requests

  - **Make sure the test suite passes** on your box, Use the provided `test_*.sh` scripts or tox.
  - Enable [Travis-Ci](http://travis-ci.org/pydata/pandas). See "Getting Travis-CI going" below.
  - Use [proper commit messages](http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html):
    - a subject line with `< 80` chars.
    - One blank line.
    - Optionally, a commit message body.
  - Please reference relevant Github issues in your commit message using `GH1234`
    or `#1234`. Either style is fine but the '#' style generates nose when your rebase your PR.
  - `doc/source/release.rst` and `doc/source/vx.y.z.txt` contain an ongoing
    changelog for each release. Add entries to these files
    as needed in a separate commit in your PR: document the fix, enhancement,
    or (unavoidable) breaking change.
  - Keep style fixes to a separate commit to make your PR more readable.
  - An informal commit message format is in effect for the project. Please try
    and adhere to it. Check `git log` for examples. Here are some common prefixes
    along with general guidelines for when to use them:
      - **ENH**: Enhancement, new functionality
      - **BUG**: Bug fix
      - **DOC**: Additions/updates to documentation
      - **TST**: Additions/updates to tests
      - **BLD**: Updates to the build process/scripts
      - **PERF**: Performance improvement
      - **CLN**: Code cleanup
  - Maintain backward-compatibility. Pandas has lots of users with lots of existing code. Don't break it.
    - If you think breakage is required clearly state why as part of the PR.
    - Be careful when changing method signatures.
    - Add deprecation warnings where needed.
  - Performance matters. Make sure your PR hasn't introduced perf regressions by using `test_perf.sh`.
  - Docstrings follow the [numpydoc](https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt) format.
  - When writing tests, use 2.6 compatible `self.assertFoo` methods. Some polyfills such as `assertRaises`
    can be found in `pandas.util.testing`.
  - Generally, pandas source files should not contain attributions. You can include a "thanks to..."
    in the release changelog. The rest is `git blame`/`git log`.
  - For extra brownie points, you can squash and reorder the commits in your PR using `git rebase -i`.
    Use your own judgment to decide what history needs to be preserved. If git frightens you, that's OK too.
  - Use `raise AssertionError` over `assert` unless you want the assertion stripped by `python -o`.
  - **Don't** merge upstream into a branch you're going to submit as a PR.
    This can create all sorts of problems. Use `git rebase` instead. This ensures
    no merge conflicts occur when your code is merged by the core team.
  - The pandas copyright policy is detailed in the pandas [LICENSE](https://github.com/pydata/pandas/blob/master/LICENSE).
  - On the subject of [PEP8](http://www.python.org/dev/peps/pep-0008/): yes.
  - On the subject of a massive PEP8-storm touching everything: not too often (once per release works).

### Notes on plotting function conventions

https://groups.google.com/forum/#!topic/pystatsmodels/biNlCvJPNNY/discussion

### Getting Travis-CI going

Instructions for getting Travis-CI installed are available [here](http://about.travis-ci.org/docs/user/getting-started/).
For those users who are new to Travis-CI and [continuous integration](https://en.wikipedia.org/wiki/Continuous_integration) in particular,
Here's a few high-level notes:
- Travis-CI is a free service (with premium account upgrades available) that integrates
  well with GitHub.
- Enabling Travis-CI on your GitHub fork of a project will cause any *new* commit
  pushed to the repo to trigger a full build+test on Travis-CI's servers.
- All the configuration for Travis-CI builds is already specified by `.travis.yml` in the repo.
  That means all you have to do is enable Travis-CI once, and then just push commits
  and you'll get full testing across py2/3 with pandas' considerable
  [test-suite](https://github.com/pydata/pandas/tree/master/pandas/tests).
- Enabling Travis-CI will attach the test results (red/green) to the Pull-Request
  page for any PR you submit. For example:

    https://github.com/pydata/pandas/pull/2532,

See the Green "Good to merge!" banner? that's it.

It's important to get travis working as PRs won't generally get merged until travis is green.

#### Steps to enable Travis-CI

- Open https://travis-ci.org/
- Select "Sign in with GitHub" (Top Navbar)
- Select \[your username\] -> "Accounts" (Top Navbar)
- Select 'Sync now' to refresh the list of repos from your GH account.
- Flip the switch for the repos you want Travis-CI enabled for.
  "pandas", obviously.
- Then, pushing a *new* commit to a certain branch on that repo
  will trigger a build/test for that branch. For example, the branch
  might be `master` or `PR1234_fix_everything__atomically`, if that's the
  name of your PR branch.

You can see the build history and current builds for your fork
at: https://travis-ci.org/(your_GH_username)/pandas.

For example, the builds for the main pandas repo can be seen at:
https://travis-ci.org/pydata/pandas.

####More developer docs

* See the [developers](http://pandas.pydata.org/developers.html) page on the
  project website for more details.
* [`pandas` wiki](https://github.com/pydata/pandas/wiki)
* [Tips and tricks](https://github.com/pydata/pandas/wiki/Tips-&-Tricks)
* [Git tips and tricks](https://github.com/pydata/pandas/wiki/Using-Git)
* [Testing advice and best practices in `pandas`](https://github.com/pydata/pandas/wiki/Testing)
