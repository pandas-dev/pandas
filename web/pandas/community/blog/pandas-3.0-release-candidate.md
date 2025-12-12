Title: pandas 3.0.0 release candidate ready for testing
Date: 2025-12-12

# pandas 3.0.0 release candidate ready for testing!

We're excited to announce the release candidate for pandas 3.0. This major
release brings significant improvements to pandas, but also features some
potentially breaking changes.

## Highlights of pandas 3.0

pandas 3.0 introduces several major enhancements:

- **Dedicated string data type by default**: String columns are now inferred as
  the new `str` dtype instead of `object`, providing better performance and type
  safety
- **Consistent copy/view behaviour with Copy-on-Write (CoW)** (a.k.a. getting
  rid of the SettingWithCopyWarning): More predictable and consistent behavior
  for all operations, with improved performance through avoiding unnecessary
  copies
- **New `pd.col` syntax**: Initial support for `pd.col()` as a simplified syntax
  for creating callables in `DataFrame.assign`
- **Enhanced deprecation policy**: A new 3-stage deprecation process to give
  downstream packages more time to adapt

You can find the complete list of changes in our
[release notes](https://pandas.pydata.org/docs/dev/whatsnew/v3.0.0.html).

## Important changes requiring code updates

As a major release, pandas 3.0 includes some breaking changes that may require
updates to your code. The two most significant changes are:

### 1. Dedicated string data type by default

Starting with pandas 3.0, string columns are automatically inferred as `str`
dtype instead of the numpy `object` (which can store any Python object).

**Example of the change:**
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

This change improves performance and type safety, but may require code updates.
For more details, see the
[String Data Type Migration Guide](https://pandas.pydata.org/docs/dev/user_guide/migration-3-strings.html).

### 2. Consistent copy/view behaviour with Copy-on-Write (CoW)

Copy-on-Write is now the default and only mode in pandas 3.0. This makes
behavior more consistent and predictable, but requires updates to certain coding
patterns:

**What you need to update:**
- **Chained assignment** will no longer work. You'll need to use `.loc` or
  `.iloc` directly on the DataFrame instead
- The `SettingWithCopyWarning` is removed (since chained assignment no longer works)
- Any indexing operation now always behaves as if it were a copy, so
  modifications won't affect the original DataFrame

**Example of the change:**
```python
# Old behavior (pandas < 3.0) - chained assignment
df[df['A'] > 0]['B'] = 1  # This might modify df (unpredictable)

# New behavior (pandas 3.0) - must use .loc
df.loc[df['A'] > 0, 'B'] = 1  # This is the correct way
```

[Copy-on-Write Migration Guide](https://pandas.pydata.org/pandas-docs/version/3.0.0/user_guide/copy_on_write.html)

## Call to Action: Test the Release Candidate

We need your help to ensure a smooth pandas 3.0 release!

Especially if you have pandas code in production or maintain a library with
pandas as a dependency, it is strongly recommended to run your test suites with
the release candidate, and report any issue to our issue tracker before the
official 3.0.0 release.

1. **Install the release candidate** and test it with your codebase
2. **Run your existing code** to identify any issues or needed updates
3. **Report any problems** you encounter on our [GitHub repository](https://github.com/pandas-dev/pandas/issues)
4. **Share your migration experiences** with the community

The more testing we get now, the smoother the final pandas 3.0 release will be
for everyone. Your feedback is crucial to making this a successful release!

### Getting the Release Candidate

You can install the pandas 3.0 release candidate from PyPI:

```bash
python -m pip install --upgrade pandas==3.0.0rc0
```

Or from conda-forge using conda/mamba:

```bash
conda install -c conda-forge/label/pandas_rc pandas==3.0.0rc0
```
