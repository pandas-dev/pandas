## Bug reports and enhancement requests

Bug reports are an important part of making *pandas* more stable. Having
a complete bug report will allow others to reproduce the bug and provide
insight into fixing. See [this stackoverflow article](https://stackoverflow.com/help/mcve)
and [this blogpost](http://matthewrocklin.com/blog/work/2018/02/28/minimal-bug-reports)
for tips on writing a good bug report.

Trying the bug-producing code out on the *master* branch is often a
worthwhile exercise to confirm the bug still exists. It is also worth
searching existing bug reports and pull requests to see if the issue has
already been reported and/or fixed.

Bug reports must:

1.  Include a short, self-contained Python snippet reproducing the
    problem. You can format the code nicely by using [GitHub Flavored
    Markdown \<http://github.github.com/github-flavored-markdown/\>]():

```python
>>> from pandas import DataFrame
>>> df = DataFrame(...)
...
```

2.  Include the full version string of *pandas* and its dependencies.
    You can use the built-in function:

```python
>>> import pandas as pd
>>> pd.show_versions()
```

3.  Explain why the current behavior is wrong/not desired and what you
    expect instead.

The issue will then show up to the *pandas* community and be open to
comments/ideas from others.
