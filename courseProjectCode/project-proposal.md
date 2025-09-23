# Project Proposal

## Project Overview

Our course project aims to build a lightweight data analysis library that
mimics essential features of the pandas ecosystem. The library will
provide tabular data structures (similar to DataFrame and Series) and
support common operations needed by scientists and engineers working
with structured data. Major functional goals include:

- **Handling missing data:** The system should represent missing values as
  NaN, NA or NaT and propagate them through computations. This
  capability simplifies data cleaning and statistical analysis by
  preventing silent errors software.com.

- **Size mutability:** Users should be able to insert or delete columns
  and rows in data structures. Dynamic resizing is central to
  interactive analysis workflows where the shape of a table evolves as
  new information becomes available raw.githubusercontent.com.

- **Automatic and explicit data alignment:** When performing
  arithmetic or merging operations, the system will align data on
  labels or allow users to opt out of alignment entirely. Proper
  alignment prevents accidental mismatches and promotes reproducible
  results raw.githubusercontent.com.

- **Flexible group-by operations:** The library should implement
  split–apply–combine patterns for aggregation, transformation, and
  filtering so that users can summarise data by categories with a
  single fluent expression raw.githubusercontent.com.

- **Robust I/O tooling:** Data structures must load from and save to
  common file formats (CSV, Excel) and efficiently persist to
  high-performance formats such as HDF5 raw.githubusercontent.com.

- **Time-series functionality:** Operations like date-range generation,
  frequency conversion, moving-window statistics and date shifting will
  be built in so that time-indexed data can be analysed without
  external libraries raw.githubusercontent.com.

In addition to these functional requirements, the project emphasises
non-functional qualities such as performance, flexibility and
expressive APIs. The goal is to provide an intuitive open-source tool
that researchers can use to analyse data without sacrificing speed or
power raw.githubusercontent.com.

---

## Key Quality Metrics

To ensure that the implementation is maintainable and testable, we will
track several quality metrics throughout the project lifecycle. The
metrics were selected based on guidance from software engineering
literature and industry best practices.

### Maintainability metrics

- **Maintainability index (MI):** Visual Studio defines an index from
  0 to 100 that summarises the ease of maintaining a piece of code.
  Higher values indicate more maintainable code, with scores above
  20 considered “good,” 10–19 “moderate” and below 10 “poor”
  learn.microsoft.com.  
  MI combines several measurements such as cyclomatic complexity,
  depth of inheritance and class coupling. Although we do not
  compute MI directly, we monitor its constituent components to track
  trends over time.

- **Cyclomatic complexity:** This measures the number of linearly
  independent paths through a program. Each decision point (e.g.,
  if, for, while) adds one to the count. Higher complexity
  indicates more potential execution paths and requires more tests to
  achieve full coverage learn.microsoft.com. Our metrics script
  approximates cyclomatic complexity by scanning for decision
  keywords, providing a reproducible indicator of structural
  complexity.

- **Comment-to-code ratio:** The number of comment lines divided by
  the number of executable lines software.com. Comments
  capture design assumptions, edge cases and rationale that are not
  obvious from code alone. A moderate ratio improves maintainability
  by aiding knowledge transfer and reducing ramp-up time for new
  contributors software.com. However, excessively high
  ratios can reflect commented-out code or verbose documentation,
  so the ratio should be interpreted in context software.com.

- **Average function length:** Smaller functions tend to perform a
  single task, are easier to understand and thus easier to modify.
  The metrics script measures the average number of code lines per
  function. Keeping this metric low encourages modular design and
  aligns with the Single Responsibility Principle.

- **Class coupling and depth of inheritance:** Although our project
  uses primarily functions and data structures, we will monitor
  coupling and inheritance depth where applicable. Visual Studio’s
  guidance notes that high class coupling and deep inheritance trees
  decrease maintainability learn.microsoft.com. We will
  minimise dependencies between modules and favour composition over
  inheritance to keep these metrics low.

### Testability metrics

- **Test coverage:** Atlassian describes code coverage as a measure
  of how much of the code base is exercised by tests and notes
  several metrics: function, statement, branch, condition and line
  coverage atlassian.com. Although a high coverage
  percentage does not guarantee good tests, it reveals which parts of
  the system remain untested and helps to prioritise additional
  testing efforts. Since we cannot run external coverage tools in
  this environment, our metrics script approximates test effort by
  reporting the ratio of lines in test files to total lines of code
  and counting the number of test functions. Increasing the
  test-to-code ratio over time should correlate with improved
  coverage.

- **Number of test cases:** We treat each test_* function in
  Python and calls to describe/it in JavaScript as individual
  test cases. Tracking the number of test cases encourages
  developers to write focused, granular tests and highlights
  subsystems that may need additional verification.

- **Complexity vs. tests:** Cyclomatic complexity informs us how
  many test cases are theoretically required to exercise all
  execution paths learn.microsoft.com. By comparing the number
  of test cases to the aggregate complexity of the code base, we can
  judge whether testing is keeping pace with growing code
  intricacy. If complexity rises faster than test counts, there may
  be untested paths that warrant attention.

---

## Using the metrics

The `metrics_collector.py` script in `courseProjectCode/Metrics/`
implements the measurements described above. Running the script
generates a JSON report containing per-file metrics and a summary.
These metrics will form the basis of our quality dashboard and guide
refactoring and testing priorities throughout the project. By
monitoring comment ratios, function lengths, complexity and test
ratios, we can make data-driven decisions to keep the code base
maintainable and to ensure that behaviour is thoroughly validated.
