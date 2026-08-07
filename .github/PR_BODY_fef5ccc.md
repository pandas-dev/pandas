## Summary
Fixes a bug in `pandas/_libs/src/parser/tokenizer.cpp` where error/warning message formatting could write through a null pointer when memory allocation fails.

This patch adds null checks after `malloc` in tokenizer error-handling paths and only calls `snprintf` when the destination pointer is non-null.

## Reproduction
Before:
- In allocation-failure paths, code could call `snprintf` with a null destination pointer.

After:
- `snprintf` is guarded by `if (ptr != NULL)` checks, so allocation-failure paths no longer attempt null-pointer writes.

## Validation
- Reviewed all touched allocation sites in tokenizer error/warning paths.
- Verified behavior is unchanged in normal (successful-allocation) paths.
- Ran relevant parser/tokenizer tests and project checks.

- [ ] closes #xxxx (Replace xxxx with the GitHub issue number)
- [ ] [Tests added and passed](https://pandas.pydata.org/pandas-docs/dev/development/contributing_codebase.html#writing-tests) if fixing a bug or adding a new feature
- [ ] All [code checks passed](https://pandas.pydata.org/pandas-docs/dev/development/contributing_codebase.html#pre-commit).
- [ ] Added [type annotations](https://pandas.pydata.org/pandas-docs/dev/development/contributing_codebase.html#type-hints) to new arguments/methods/functions.
- [ ] Added an entry in the latest `doc/source/whatsnew/vX.X.X.rst` file if fixing a bug or adding a new feature.
- [x] I have reviewed and followed all the [contribution guidelines](https://pandas.pydata.org/docs/dev/development/contributing.html)
- [x] If I used AI to develop this pull request, I prompted it to follow `AGENTS.md`.
