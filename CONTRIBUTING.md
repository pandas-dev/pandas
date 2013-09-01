###Guidelines

All contributions, bug reports, bug fixes, documentation improvements,
enhancements and ideas are welcome.

The Github "issues" tab contains some issues labels "Good as first PR", these are
tasks which do not require deep knowledge of the package. Look those up if you're
looking for a quick way to help out.

Please try and follow these guidelines, as this makes it easier for us to accept
your contribution or address the issue you're having.

- When submitting a bug report:
  - Please include a short, self-contained python snippet reproducing the problem.
  You can have the code formatted nicely by using [GitHub Flavored Markdown](http://github.github.com/github-flavored-markdown/) :

    ```
    ```python

    print("I ♥ pandas!")

    ``'
    ```

  - Specify the pandas (and numpy) version used. (you can check `pandas.__version__`).
  - Explain what the expected behavior was, and what you saw instead.
  - If the issue seems to involve some of pandas' dependencies such as matplotlib
    or PyTables, you should include (the relavent parts of) the output of
    [ci/print_versions.py](https://github.com/pydata/pandas/blob/master/ci/print_versions.py)

- When submitting a Pull Request
  - **Make sure the test suite passes**., and that means on python3 as well.
    You can use "test_fast.sh", or tox locally and/or enable Travis-CI on your fork.
    See the "Getting Travis-CI going" below.
  - If you are changing any code, you need to enable Travis-CI on your fork,
    to make it easier for the team to see that the PR does indeed pass all the tests.
  - Back-compatibility **really** matters. Pandas already has a large user-base and
    a lot of existing user code. Don't break old code if you can avoid it
    Explain the need if there is one in the PR.
    Changes to method signatures should be made in a way which doesn't break existing
    code, for example you should beware of changes to ordering and naming of keyword
    arguments. Add deprecation warnings when needed.
  - Performance matters. You can use the included "test_perf.sh"
    script to make sure your PR does not introduce any performance regressions
    in the library.
  - docstrings follow the [numpydoc](https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt) format.
  - **Don't** merge upstream into a branch you're going to submit as a PR,
    This can create all sorts of problems. Use "git rebase" instead. This ensures
    no merge conflicts occur when you're code is merged by the core team.
  - An informal commit message format is in effect for the project, please try
    and adhere to it. View "git log" for examples. Here are some common prefixes
    along with general guidelines for when to use them:
      - ENH: Enhancement, new functionality
      - BUG: Bug fix
      - DOC: Additions/updates to documentation
      - TST: Additions/updates to tests
      - BLD: Updates to the build process/scripts
      - PERF: Performance improvement
      - CLN: Code cleanup
  - Commit messages should have subject line <80 chars, followed by one blank line,
    and finally a commit message body if there's a need for one.
  - Please reference the GH issue number in your commit message using GH1234
    or #1234, either style is fine.
  - Use "raise AssertionError" rather then plain `assert` in library code (using assert is fine
    for test code). python -o strips assertions. better safe then sorry.
  - When writing tests, don't use "new" assertion methods added to the unittest module
    in 2.7 since pandas currently supports 2.6. The most common pitfall is:

    with self.assertRaises(ValueError):
         foo

    which fails on python 2.6. You need to use `assertRaises` from
    `pandas.util.testing` instead (or use `self.assertRaises(TheException,func,args)`).

  - doc/source/release.rst and doc/source/vx.y.z.txt contain an on-going
    changelog for each release as it is worked on. Add entries to these files
    as needed in a separate commit in your PR, documenting the fix, enhancement
    or (unavoidable) breaking change.
  - For extra brownie points, use "git rebase -i" to squash and reorder
    commits in your PR so that the history makes the most sense. Use your own
    judgment to decide what history needs to be preserved.
  - Pandas source code should not (with some exceptions, such as 3rd party licensed code),
    generally speaking, include an "Authors:" list or attribution to individuals in source code.
    The RELEASE.rst details changes and enhancements to the code over time,
    a "thanks goes to @JohnSmith." as part of the appropriate entry is a suitable way to acknowledge
    contributions, the rest is git blame/log.
    Feel free to ask the commiter who merges your code to include such an entry
    or include it directly yourself as part of the PR if you'd like to. We're always glad to have
    new contributors join us from the ever-growing pandas community.
    You may also be interested in the copyright policy as detailed in the pandas [LICENSE](https://github.com/pydata/pandas/blob/master/LICENSE).
  - On the subject of [PEP8](http://www.python.org/dev/peps/pep-0008/): yes.
  - On the subject of massive PEP8 fix PRs touching everything, please consider the following:
    - They create merge conflicts for people working in their own fork.
    - They makes git blame less effective.
    - Different tools / people achieve PEP8 in different styles. This can create
      "style wars" and churn that produces little real benefit.
    - If your code changes are intermixed with style fixes, they are harder to review
      before merging. Keep style fixes in separate commits.
    - it's fine to clean-up a little around an area you just worked on.
    - Generally its a BAD idea to PEP8 on documentation.

    Having said that, if you still feel a PEP8 storm is in order, go for it.

### Notes on plotting functions convention

https://groups.google.com/forum/#!topic/pystatsmodels/biNlCvJPNNY/discussion

###Getting Travis-CI going

Instructions for getting Travis-CI installed are available [here](http://about.travis-ci.org/docs/user/getting-started/). For those users who are new to travis-ci and continuous-integration in particular,
Here's a few high-level notes:
- Travis-CI is a free service (with premium account available), that integrates
well with Github.
- Enabling Travis-CI on your github fork of a project will cause any *new* commit
pushed to the repo to trigger a full build+test on Travis-CI's servers.
- All the configuration for travis's builds is already specified by .travis.yml in the repo,
That means all you have to do is enable Travis-CI once, and then just push commits
and you'll get full testing across py2/3 with pandas's considerable test-suite.
- Enabling travis-CI will attach the test results (red/green) to the Pull-Request
page for any PR you submit. For example:

    https://github.com/pydata/pandas/pull/2532,

See the Green "Good to merge!" banner? that's it.

This is especially important for new contributors, as members of the pandas dev team
like to know the test suite passes before considering it for merging.
Even regular contributors who test religiously on their local box (using tox
for example) often rely on a PR+travis=green to make double sure everything
works ok on another system, as occasionally, it doesn't.

####Steps to enable Travis-CI

- go to https://travis-ci.org/
- "Sign in with Github", on top panel.
- \[your username\]/Account, on top-panel.
- 'sync now' to refresh the list of repos on your GH account.
- flip the switch on the repos you want Travis-CI enabled for,
"pandas" obviously.
- Then, pushing a *new* commit to a certain branch on that repo
will trigger a build/test for that branch, for example the branch
might be "master" or "PR1234_fix_all_the_things", if that's the
name of your PR branch.

You can see the build history and current builds for your fork
on: https://travis-ci.org/(your_GH_username)/pandas.

For example, the builds for the main pandas repo can be seen at:
https://travis-ci.org/pydata/pandas.

####More developer docs

* See the [developers](http://pandas.pydata.org/developers.html) page on the
  project website for more details.
* [`pandas` wiki](https://github.com/pydata/pandas/wiki)
* [Tips and tricks](https://github.com/pydata/pandas/wiki/Tips-&-Tricks)
* [Git tips and tricks](https://github.com/pydata/pandas/wiki/Using-Git)
* [Testing advice and best practices in `pandas`](https://github.com/pydata/pandas/wiki/Testing)
