# Security Policy

This document describes the security model of pandas and how to report
vulnerabilities.

## Security Model

pandas is a data analysis library. It reads, transforms, and writes data on
behalf of the user, with the permissions of the Python process that runs it.
pandas does not provide a sandbox: code that calls pandas is trusted, and
arguments passed to pandas functions (file paths, URLs, expressions, SQL
queries, options) are trusted.

Functions that read data are intended to be usable with data from an
untrusted source, except for the areas listed under
[Areas Outside the Security Model](#areas-outside-the-security-model). For
functions outside those areas, we consider an issue a **security
vulnerability** if an attacker who controls the data could exploit it to:

- Execute arbitrary code (Remote Code Execution); or
- Read sensitive information from process memory (Information Disclosure).

Other unexpected behavior caused by malformed or malicious data is generally
considered a **bug**, not a security vulnerability. This includes:

- Crashes (including segmentation faults), exceptions, or incorrect results.
- Excessive memory or CPU use, very long running time, or infinite loops.

If you are unsure whether an issue is a security vulnerability, report it
privately as described in [Reporting a Vulnerability](#reporting-a-vulnerability).

## Areas Outside the Security Model

We occasionally receive vulnerability reports for the areas below. pandas does
not provide a security boundary in these areas, so we are unlikely to consider
such reports security vulnerabilities. If you are unsure whether an issue
in one of these areas is a security vulnerability, report it privately and we
will discuss it.

### Pickle

`read_pickle` and `DataFrame.to_pickle` use Python's
[pickle](https://docs.python.org/3/library/pickle.html) module, which can
execute arbitrary Python code when loading data. `read_hdf` and `HDFStore` also
use pickle to store columns with `object` dtype. pandas does not provide any
security on top of pickle. Only load pickled data that you trust.

### Expression evaluation

`eval`, `DataFrame.eval`, and `DataFrame.query` evaluate string expressions in
the context of a DataFrame using numexpr or Python. pandas does not provide any
security on top of these engines. Do not pass expressions from untrusted
sources to these methods.

### XSLT stylesheets

`read_xml` and `DataFrame.to_xml` accept a `stylesheet` argument that runs an
XSLT program through lxml. pandas does not provide any security on top of lxml.
Do not use stylesheets from untrusted sources.

### HTML output

`DataFrame.to_html` escapes HTML characters by default (`escape=True`), and
`Styler` can escape HTML when `escape="html"` is set (it does not by default).
These options exist for convenience and are not a security boundary. pandas
does not guarantee that the HTML it generates is safe to serve to a browser
when the data or the formatting options come from an untrusted source.

### Injection into other formats

Several methods read or write formats where values can be interpreted as code
or commands, such as SQL (`read_sql`) and spreadsheet formulas (`to_csv`,
`to_excel`). pandas does not escape or sanitize values for these formats.
Use the parameter binding or escaping options of the underlying library
(for example, `params` in `read_sql`) for untrusted values.

### File paths and URLs

pandas reads from and writes to any file path or URL it is given, including
remote locations through fsspec. pandas does not restrict which paths or URLs
can be accessed. Applications that pass untrusted paths or URLs to pandas must
validate them first.

### Clipboard

`read_clipboard` and `DataFrame.to_clipboard` run external programs (such as
`xclip`, `xsel`, `wl-copy`, or `pbcopy`) found on the `PATH`.

### Third-party libraries

pandas uses many third-party libraries, such as NumPy, PyArrow, lxml,
openpyxl, and SQLAlchemy. Please report vulnerabilities in those libraries to
their maintainers. A report is in scope for pandas only if pandas itself uses
the library in an unsafe way.

## Supported Versions

Security fixes are made in the latest release of pandas. Fixes are not
backported to older minor versions.

## Reporting a Bug

We take all bugs seriously and welcome help fixing them. If you find a bug that
clearly does not meet the criteria for a security vulnerability above, please
report it in the [public issue tracker](https://github.com/pandas-dev/pandas/issues).

## Reporting a Vulnerability

**Do not report security vulnerabilities in a public issue, pull request, or
discussion.**

Report security vulnerabilities privately through
[GitHub private vulnerability reporting](https://github.com/pandas-dev/pandas/security/advisories/new).
This allows the maintainers to investigate and release a fix before the
vulnerability is made public.

Please include in your report:

- A clear description of the issue and a minimal reproducer.
- The affected pandas version, and the versions of related libraries
  (for example, NumPy or PyArrow) and the operating system.
- The potential impact, including how an attacker could exploit it.

## Security Advisories

When we confirm a security vulnerability, we publish a
[GitHub security advisory](https://github.com/pandas-dev/pandas/security/advisories)
and request a CVE as appropriate. Issues that are bugs under the security
model above, such as crashes that are not exploitable, are fixed as ordinary
bugs and do not receive an advisory.
