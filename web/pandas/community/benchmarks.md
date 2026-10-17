# Benchmarks

Benchmarks are tests to measure the performance of pandas. There are two different
kinds of benchmarks relevant to pandas:

* Internal pandas benchmarks to measure speed and memory usage over time
* Community benchmarks comparing the speed or memory usage of different tools at
  doing the same job

## pandas benchmarks

pandas benchmarks are implemented in the [asv_bench](https://github.com/pandas-dev/pandas/tree/main/asv_bench)
directory of our repository, using the [airspeed velocity](https://asv.readthedocs.io/en/latest/)
(asv for short) framework. See
[running the performance test suite]({{ base_url }}docs/development/contributing_codebase.html#running-the-performance-test-suite)
in the contributing guide for how to run them.

## Automated benchmark runner

The [asv-runner](https://github.com/pandas-dev/asv-runner/) repository automatically runs the pandas asv benchmark suite
for every (or almost every) commit to the `main` branch. It is run on GitHub actions.
See the linked repository for more details. The results are available at:

[https://pandas-dev.github.io/asv-runner/](https://pandas-dev.github.io/asv-runner/)

## Community benchmarks

The main benchmarks comparing dataframe tools that include pandas are:

- [DuckDB (former H2O.ai) benchmarks](https://duckdblabs.github.io/db-benchmark/)
- [TPCH benchmarks](https://pola.rs/posts/benchmarks/)
