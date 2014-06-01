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

  - Include the full version string of pandas and it's dependencies. In recent (>0.12) versions
    of pandas you can use a built in function:

    ```python
    >>> from pandas.util.print_versions import show_versions
    >>> show_versions()
    ```

    and in 0.13.1 onwards:
    ```python
    >>> pd.show_versions()
    ```

  - Explain what the expected behavior was, and what you saw instead.

#### Pull Requests

  - **Make sure the test suite passes** on your box, Use the provided `test_*.sh` scripts or tox.
  - Use [proper commit messages](http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html):
    - a subject line with `< 80` chars.
    - One blank line.
    - Optionally, a commit message body.
  - Please reference relevant Github issues in your commit message using `GH1234`
    or `#1234`. Either style is fine but the '#' style generates noise when your rebase your PR.
  - `doc/source/vx.y.z.txt` contains an ongoing
    changelog for each release. Add an entry to this file
    as needed in your PR: document the fix, enhancement,
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
  - Write tests.
  - When writing tests, use 2.6 compatible `self.assertFoo` methods. Some polyfills such as `assertRaises`
    can be found in `pandas.util.testing`.
  - Do not attach doctrings to tests. Make the test itself readable and use comments if needed.
  - Generally, pandas source files should not contain attributions. You can include a "thanks to..."
    in the release changelog. The rest is `git blame`/`git log`.
  - When you start working on a PR, start by creating a new branch pointing at the latest
    commit on github master.
  - **Do not** merge upstream into a branch you're going to submit as a PR.
    Use `git rebase` against the current github master.
  - For extra brownie points, you can squash and reorder the commits in your PR using `git rebase -i`.
    Use your own judgment to decide what history needs to be preserved. If git frightens you, that's OK too.
  - Use `raise AssertionError` over `assert` unless you want the assertion stripped by `python -o`.
  - The pandas copyright policy is detailed in the pandas [LICENSE](https://github.com/pydata/pandas/blob/master/LICENSE).
  - On the subject of [PEP8](http://www.python.org/dev/peps/pep-0008/): yes.
  - We've written a tool to check that your commits are PEP8 great,
    [`pip install pep8radius`](https://github.com/hayd/pep8radius). Look at PEP8 fixes in your branch
    vs master with `pep8radius master --diff` and make these changes with
    `pep8radius master --diff --in-place`.
  - On the subject of a massive PEP8-storm touching everything: not too often (once per release works).

### Notes on plotting function conventions

https://groups.google.com/forum/#!topic/pystatsmodels/biNlCvJPNNY/discussion

####More developer docs

* See the [developers](http://pandas.pydata.org/developers.html) page on the
  project website for more details.
* [`pandas` wiki](https://github.com/pydata/pandas/wiki)
* [Tips and tricks](https://github.com/pydata/pandas/wiki/Tips-&-Tricks)
* [Git tips and tricks](https://github.com/pydata/pandas/wiki/Using-Git)
* [Testing advice and best practices in `pandas`](https://github.com/pydata/pandas/wiki/Testing)
