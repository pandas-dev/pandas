# PDEP-4: Consistent datetime parsing

- Created: 18 September 2022
- Status: Implemented
- Discussion: [#48621](https://github.com/pandas-dev/pandas/pull/48621)
- Author: [Marco Gorelli](https://github.com/MarcoGorelli)
- Revision: 2

[TOC]

## Abstract

The suggestion is that:

- ``to_datetime`` becomes strict and uses the same datetime format to parse all elements in its input.
  The format will either be inferred from the first non-NaN element (if `format` is not provided by the user), or from
  `format`;
- ``infer_datetime_format`` be deprecated (as a strict version of it will become the default);
- an easy workaround for non-strict parsing be clearly documented.

## Motivation and Scope

Pandas date parsing is very flexible, but arguably too much so - see
https://github.com/pandas-dev/pandas/issues/12585 and linked issues for how
much confusion this causes. Pandas can swap format midway, and though this
is documented, it regularly breaks users' expectations.

Simple example:
```ipython
In [1]: pd.to_datetime(['12-01-2000 00:00:00', '13-01-2000 00:00:00'])
Out[1]: DatetimeIndex(['2000-12-01', '2000-01-13'], dtype='datetime64[ns]', freq=None)
```
The user was almost certainly intending the data to be read as "12th of January, 13th of January".
However, it's read as "1st of December, 13th of January". No warning or error is thrown.

Currently, the only way to ensure consistent parsing is by explicitly passing
``format=``. The argument ``infer_datetime_format``
isn't strict, can be called together with ``format``, and can still break users' expectations:

```ipython
In [2]: pd.to_datetime(['12-01-2000 00:00:00', '13-01-2000 00:00:00'], infer_datetime_format=True)
Out[2]: DatetimeIndex(['2000-12-01', '2000-01-13'], dtype='datetime64[ns]', freq=None)
```

## Detailed Description

Concretely, the suggestion is:

- if no ``format`` is specified, ``pandas`` will guess the format from the first non-NaN row
  and parse the rest of the input according to that format. Errors will be handled
  according to the ``errors`` argument - there will be no silent switching of format;
- ``infer_datetime_format`` will be deprecated;
- ``dayfirst`` and ``yearfirst`` will continue working as they currently do;
- if the format cannot be guessed from the first non-NaN row, a ``UserWarning`` will be thrown,
  encouraging users to explicitly pass in a format.
  Note that this should only happen for invalid inputs such as `'a'`
  (which would later throw a ``ParserError`` anyway), or inputs such as ``'00:12:13'``,
  which would currently get converted to ``''2022-09-18 00:12:13''``.

If a user has dates in a mixed format, they can still use flexible parsing and accept
the risks that poses, e.g.:
```ipython
In [3]: pd.to_datetime(['12-01-2000 00:00:00', '13-01-2000 00:00:00'], format='mixed')
Out[3]: DatetimeIndex(['2000-12-01', '2000-01-13'], dtype='datetime64[ns]', freq=None)
```
or, if their dates are all ISO8601,
```ipython
In [4]: pd.to_datetime(['2020-01-01', '2020-01-01 03:00'], format='ISO8601')
Out[4]: DatetimeIndex(['2020-01-01 00:00:00', '2020-01-01 03:00:00'], dtype='datetime64[ns]', freq=None)
```

## Usage and Impact

My expectation is that the impact would be a net-positive:

- potentially severe bugs in people's code will be caught early;
- users who actually want mixed formats can still parse them, but now they'd be forced to be
  very explicit about it;
- the codebase would be noticeably simplified.

As far as I can tell, there is no chance of _introducing_ bugs.

## Implementation

The whatsnew notes read

> In the next major version release, 2.0, several larger API changes are being considered without a formal deprecation.

I'd suggest making this change as part of the above, because:

- it would only help prevent bugs, not introduce any;
- given the severity of bugs that can result from the current behaviour, waiting another 2 years until pandas 3.0.0
  would potentially cause a lot of damage.

Note that this wouldn't mean getting rid of ``dateutil.parser``, as that would still be used within ``guess_datetime_format``. With this proposal, however, subsequent rows would be parsed with the guessed format rather than repeatedly calling ``dateutil.parser`` and risk having it silently switch format

Finally, the function ``from pandas._libs.tslibs.parsing import guess_datetime_format`` would be made public, under ``pandas.tools``.

## Out of scope

We could make ``guess_datetime_format`` smarter by using a random sample of elements to infer the format.

### PDEP History

- 18 September 2022: Initial draft
- 25 January 2023: Amended to mention ``format='ISO8601'`` and ``format='mixed'`` options
