# Benchmarks

Benchmarks are tests to measure the performance of pandas. There are two different
kinds of benchmarks relevant to pandas:

* Internal pandas benchmarks to measure speed and memory usage over time
* Community benchmarks comparing the speed or memory usage of different tools at
  doing the same job

## pandas benchmarks

pandas benchmarks are implemented in the [asv_bench](https://github.com/pandas-dev/pandas/tree/main/asv_bench)
directory of our repository. The benchmarks are implemented for the
[airspeed velocity](https://asv.readthedocs.io/en/latest/) (asv for short) framework.

The benchmarks can be run locally by any pandas developer. This can be done
with the `asv run` command, and it can be useful to detect if local changes have
an impact in performance, by running the benchmarks before and after the changes.
More information on running the performance test suite is found
[here](https://pandas.pydata.org/docs/dev/development/contributing_codebase.html#running-the-performance-test-suite).

Note that benchmarks are not deterministic, and running in different hardware or
running in the same hardware with different levels of stress have a big impact in
the result. Even running the benchmarks with identical hardware and almost identical
conditions can produce significant differences when running the same exact code.

## Automated benchmark runner

The [asv-runner](https://github.com/pandas-dev/asv-runner/) repository automatically runs the pandas asv benchmark suite
for every (or almost every) commit to the `main` branch. It is run on GitHub actions.
See the linked repository for more details. The results are available at:

https://pandas-dev.github.io/asv-runner/

## Community benchmarks

The main benchmarks comparing dataframe tools that include pandas are:

- [DuckDB (former H2O.ai) benchmarks](https://duckdblabs.github.io/db-benchmark/)
- [TPCH benchmarks](https://pola.rs/posts/benchmarks/)
