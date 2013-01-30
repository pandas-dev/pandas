Guidelines
---

All contributions, bug reports, bug fixes, documentation improvments,
enhancements and ideas are welcome. Please try and follow these guidelines
as best you can, this makes it easier for us to accept your contribution or
address the issue you're having.

- When submitting a bug report:
  - Please include a short, self-contained python snippet.
  - Specify the pandas version used. (you can check `pandas.__version__`).
  - Explain what the expected behavior was, and what you saw instead.
- When submitting a Pull Request
  - **Make sure the test suite passes**., and that means on python3 as well.
    You can use "test_fast.sh", or tox locally and/or enable Travis-CI on you're fork.
  - We suggest you enable Travis-CI on your fork, to make it easier for the team
     to see that the PR does indeed pass all the tests.
  - Pandas cares about performance, you can use the included "test_perf.sh"
    script to make sure your PR does not introduce any performance regressions
    in the library.
  - **Don't** merge upstream into a branch you're going to submit as a PR,
    This can create all sorts of problems. Use "git rebase" instead. This ensures
    no merge conflicts occur when you're code is merged by the core team.
  - An informal commit message format is in use by the project, please try
    and adhere to it. Use a "ENH: ", "TST:", "BUG:", "DOC:", etc' prefix in
    your commit title. Check the output of "git log" for examples.
  - doc/source/vx.y.z.txt contains an on-going change-log for the version under
    development, look at previous versions, and add an entry as needed in
    a separate commit in your PR.
  - For extra brownie points, use "git rebase -i" to squash and reorder
    commits in your PR so that the history makes the most sense. Use your own
    judgment to decide what history needs to be preserved.

Please see [Developers](http://pandas.pydata.org/developers.html) page on
the project website for more details.
