Title: pandas 3.0 released!
Date: 2026-01-21

# pandas 3.0 released!

We're excited to announce the release of pandas 3.0.0. This major
long-awaited release brings significant improvements to pandas, but also
features some potentially breaking changes.

## Highlights of pandas 3.0

pandas 3.0 introduces several major enhancements:

- **Dedicated string data type by default**: string columns are now inferred as
  the new `str` dtype instead of `object`, providing better performance and type
  safety
- **Consistent copy/view behaviour with Copy-on-Write (CoW)** (a.k.a. getting
  rid of the SettingWithCopyWarning): more predictable and consistent behavior
  for all operations, with improved performance through avoiding unnecessary
  copies
- **New default resolution for datetime-like data**: no longer defaulting to
  nanoseconds, but generally microseconds (or the resolution of the input), when
  when constructing datetime or timedelta data (avoiding out-of-bounds errors
  for dates with a year before 1678 or after 2262)
- **New `pd.col` syntax**: initial support for `pd.col()` as a simplified syntax
  for creating callables in `DataFrame.assign`

Further, pandas 3.0 includes a lot of other improvements and bug fixes. You can
find the complete list of changes in our
[release notes](https://pandas.pydata.org/docs/dev/whatsnew/v3.0.0.html).

## Upgrading to pandas 3.0

As a major release, pandas 3.0 includes some breaking changes that may require
updates to your code. The two most significant changes are:

### 1. Dedicated string data type by default

Starting with pandas 3.0, string columns are automatically inferred as `str`
dtype instead of the numpy `object` (which can store any Python object).

**Example:**
```python
# Old behavior (pandas < 3.0)
>>> ser = pd.Series(["a", "b"])
>>> ser
0    a
1    b
dtype: object  # <-- numpy object dtype

# New behavior (pandas 3.0)
>>> ser = pd.Series(["a", "b"])
>>> ser.dtype
>>> ser
0    a
1    b
dtype: str  # <-- new string dtype
```

This change improves performance and type safety, but may require code updates,
especially for library code that currently looks for "object" dtype when
expecting string data.

For more details, see the
[migration guide for the new string data type](https://pandas.pydata.org/docs/dev/user_guide/migration-3-strings.html).

### 2. Consistent copy/view behaviour with Copy-on-Write (CoW)

Copy-on-Write is now the default and only mode in pandas 3.0. This makes
behavior more consistent and predictable, but requires updates to certain coding
patterns.

The most impactfull change is that **chained assignment will no longer work**.
As a result, the `SettingWithCopyWarning` is also removed (since there is no
longer ambiguity whether it would work or not), and defensive `.copy()` calls
to silence the warning are no longer needed.

**Example:**
```python
# Old behavior (pandas < 3.0) - chained assignment
df["foo"][df["bar"] > 5] =   # This might modify df (unpredictable)

# New behavior (pandas 3.0) - must do the modification in one step (e.g. with .loc)
df.loc[df["bar"] > 5, "foo"] = 100
```

In general, any result of an indexing operation or method now always behaves as
if it were a copy, so modifications of the result won't affect the original
DataFrame.

For more details, see the
[Copy-on-Write migration guide](https://pandas.pydata.org/docs/dev/user_guide/copy_on_write.html#migrating-to-copy-on-write).


## Call to Action: test the Release Candidate

We need your help to ensure a smooth pandas 3.0 release!

Especially if you have pandas code in production or maintain a library with
pandas as a dependency, it is strongly recommended to run your test suites with
the release candidate, and report any issue to our issue tracker before the
official 3.0.0 release.

How can you best test the release candidate?

1. **First update to the latest released pandas 2.3** (if you are not already
   running that version) and test it with your codebase. It is recommended to
   resolve any deprecation warning before upgrading to pandas 3.0.
2. Optionally, you can already enable the new string dtype and Copy-on-Write
   mode using pandas 2.3 (`pd.options.future.infer_string = True` and
   `pd.options.mode.copy_on_write = True`).
3. **Install the release candidate** (see below) and test it with your codebase
4. **Run your existing code** to identify any issues or needed updates
5. **Report any problems** you encounter on our [GitHub repository issue tracker](https://github.com/pandas-dev/pandas/issues)

The more testing we get now, the smoother the final pandas 3.0 release will be
for everyone. Your feedback is crucial for making this a successful release!

### Getting pandas 3.0

You can install the latest pandas 3.0 release from PyPI:

```bash
python -m pip install --upgrade pandas==3.0.*
```

Or from conda-forge using conda/mamba:

```bash
conda install -c conda-forge pandas=3.0
```
