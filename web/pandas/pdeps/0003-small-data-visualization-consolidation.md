# PDEP-3: Small Data Visualization Consolidation

- Created: 13 August 2022
- Status: Draft
- Discussion: [#XXXXX](https://github.com/pandas-dev/pandas/pull/XXXXX)
- Author: [Hamish Darbyshire](https://github.com/attack68)
- Revision: 1

## Abstract

Developing consistent and flexible output methods for small data,
which are well documented and decoupled from `DataFrame` methods.

## Motivation and Scope

Large data is not in scope for this PDEP. Large Data is classed as that which
would not be expected to print, or be wire-transferred for example, with LaTeX,
HTML or JSON.

An important part of data analysis is data exploration and data visualization.
`Styler` already exists as a part of Pandas, and whose original purpose was
for conditional formatting for data cells in HTML. A rough timeline of
activity of Styler development has been:
  - 2015: released as conditional HTML formatting
  - 2015-2020: expanded CSS functionality and added built in conditional
    formatting functions.
  - 2020: performance enhancements,
  - 2021: expanded built in methods, expanded formatting functions,
    performance enhancements, added LaTeX as an alternative to HTML,
  - 2022: added concat to combine the output of multiple Stylers, documentation
    and consistency improvements.

With the feature enhancements in 2021/22 coupled with performance improvements
it is appropriate to re-target the Styler away from an HTML conditional
formatting only tool.

The overall design philosophy for Styler is proposed as:

  - i) Function as a decoupled display formatting tool.
  - ii) Provide output methods for small data within a consistent API framework.
  - iii) Provide enhanced formatting techniques such as hiding, color, font and border control to output methods that support it, and document these.
  - iv) Should be the only displaying rendering method for HTML and LaTeX.

## Detailed Description

*Styler Philosophy - i) Decoupled display formatter:*

This involves creating output without altering, or the need to pre-alter any of
the in memory data of the DataFrame.
This serves the purpose of maintaining high performance DataFrame indexing
and data features without penalising display.

```python
>>> df = DataFrame({"data": range(100)},
                   columns=date_range("2022-07-01", periods=100))
>>> df.loc["2022-08-10": "2022-08-14"].style.format_index(lambda v: v.strftime("%a"))
           data
Wed   -0.245266
Thu    0.123456
Fri    1.123456
Sat   -1.124578
Sun    0.999888
```

It also allows concatenation of display without the need to concatenate DataFrames
and avoids mixing column dtypes.

```python
>>> df = DataFrame({"norm": np.random.randn(10000),
                "poisson": np.random.poisson(1, 10000)})
>>> df_summary = df.agg(["mean"])
>>> pd.options.styler.render.max_rows=5
>>> df.style.concat(df_summary.style)
           norm    poisson
0     -0.111111          0
1      0.111111          1
2      1.234567          0
3     -2.666777          1
4      0.443355          2
...         ...        ...
mean   0.001122   0.999999
```

This design choice allows it to be separated from core functions,
permitting a forked component in the future if ever necessary.

*Styler Philosophy - ii) Output methods for small data*

"Within a consistent framework" means that the method chaining construct
should be applicable to all output methods, for example the following should
all experience similar behaviour (as well as others):

  - `styler.hide(*args).format(*args).to_latex()`
  - `styler.hide(*args).format(*args).to_html()`
  - `styler.hide(*args).format(*args).to_json()`

Any exceptions to this rule are documented and transparent (for example Excel
cannot currently implement either `hide` or `format`. This allows for generalised
styler construction with the output method determinable at render time.

The full list of proposed output methods from Styler is
- `to_html`, (implemented with jinja2 templates: fully functional)
- `to_latex`, (implemented with jinja2 templates: fully functional)
- `to_json`, (not implemented)
- `to_string`, (implemented with jinja2 templates: basic initial commit version,
   needs extension for appropriate console printing)
- `to_csv`, (not implemented, albeit indirectly available via to_string)
- `to_excel`, (available via augmenting DataFrame.to_excel: this is not decoupled)

*Styler Philosophy  -iii) Enhanced formatting and documentation*

For the most visual of outputs, HTML and LaTeX the suite of functionality
predominantly exists. For excel progress has been towards unification
but still some features are incompatible.

Documentation development is an important aspect here to unify all
methods and give user examples.

*Styler Philosophy  -iv) Only renderer for HTML and LaTeX*

DataFrame LatexFormatter and HTMLFormatter exist for legacy implementations.
They have comparatively fewer features, are less performant by some metrics,
and in the case of HTML contain deprecated HTML, and potentially non-CSP valid
output (e.g. inline styles). The proposal is to keep `DataFrame.to_html` and
`DataFrame.to_latex` methods but redesign their arguments signature with a
view to being more consistent with the arguments of Styler, and create output
via the Styler implementation, thus with a requirement for `jinja2`.

## Usage and Impact

It is expected that these tools could feasibly be used by any user.
Small data users may have requirements to use the tools directly,
whilst big data users will often create summary tables to explore and
examine the data where this would otherwise be useful.

Providing consistent functionality across outputs and well
documented formatting features will add to the overall appeal
of the pandas package and promote its longevity, as a "single,
fully featured package".

## Implementation

A number of release notes, 1.1, 1.2, 1.3, 1.4, have already documented
development towards these objectives,

The required implementation is to:

  - advance the outstanding output methods that are not yet implemented,
    or partly implemented, or are partly conforming to the philosophy
    (e.g. Styler.to_excel)
  - synchronise the mentioned DataFrame output methods to utilise the
    Styler implementations, and alter any keyword arguments.
  - revise and restructure the documentation to present Styler as a
    holistic output formatter with superior examples of relevant formats.

### PDEP History

- 13 August 2022: Initial draft
