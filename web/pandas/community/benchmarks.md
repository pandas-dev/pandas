# Benchmarks

Benchmarks are tests to measure the performance of pandas. There are two different
kinds of benchmarks relevant to pandas:

* Internal pandas benchmarks to measure speed and memory usage over time
* Community benchmarks comparing the speed or memory usage of different tools at
  doing the same job

## pandas benchmarks

pandas benchmarks are implemented in the [asv_bench](https://github.com/pandas-dev/pandas/tree/main/asv_bench)
directory of our repository. The benchmarks are implemented for the
[airspeed velocity](https://asv.readthedocs.io/en/v0.6.1/) (asv for short) framework.

The benchmarks can be run locally by any pandas developer. This can be done
with the `asv run` command, and it can be useful to detect if local changes have
an impact in performance, by running the benchmarks before and after the changes.
More information on running the performance test suite is found
[here](https://pandas.pydata.org/docs/dev/development/contributing_codebase.html#running-the-performance-test-suite).

Note that benchmarks are not deterministic, and running in different hardware or
running in the same hardware with different levels of stress have a big impact in
the result. Even running the benchmarks with identical hardware and almost identical
conditions produces significant differences when running the same exact code.

## pandas benchmarks servers

We currently have two physical servers running the benchmarks of pandas for every
(or almost every) commit to the `main` branch. The servers run independently from
each other. The original server has been running for a long time, and it is physically
located with one of the pandas maintainers. The newer server is in a datacenter
kindly sponsored by [OVHCloud](https://www.ovhcloud.com/). More information about
pandas sponsors, and how your company can support the development of pandas is
available at the [pandas sponsors]({{ base_url }}about/sponsors.html) page.

Results of the benchmarks are available at:

- Original server: [asv](https://asv-runner.github.io/asv-collection/pandas/)
- OVH server: [asv](https://pandas.pydata.org/benchmarks/asv/) (benchmarks results can
  also be visualized in this [Conbench PoC](http://57.128.112.95:5000/)

### Original server configuration

The machine can be configured with the Ansible playbook in
[tomaugspurger/asv-runner](https://github.com/tomaugspurger/asv-runner).
The results are published to another GitHub repository,
[tomaugspurger/asv-collection](https://github.com/tomaugspurger/asv-collection).

The benchmarks are scheduled by [Airflow](https://airflow.apache.org/).
It has a dashboard for viewing and debugging the results.
You’ll need to setup an SSH tunnel to view them:

```
ssh -L 8080:localhost:8080 pandas@panda.likescandy.com
```

### OVH server configuration

The server used to run the benchmarks has been configured to reduce system
noise and maximize the stability of the benchmarks times.

The details on how the server is configured can be found in the
[pandas-benchmarks repository](https://github.com/pandas-dev/pandas-benchmarks).
There is a quick summary here:

- CPU isolation: Avoid user space tasks to execute in the same CPU as benchmarks, possibly interrupting them during the execution (include all virtual CPUs using a physical core)
- NoHZ: Stop the kernel tick that enables context switching in the isolated CPU
- IRQ affinity: Ban benchmarks CPU to avoid many (but not all) kernel interruption in the isolated CPU
- TurboBoost: Disable CPU scaling based on high CPU demand
- P-States: Use "performance" governor to disable P-States and CPU frequency changes based on them
- C-States: Set C-State to 0 and disable changes to avoid slower CPU after system inactivity

## Community benchmarks

The main benchmarks comparing dataframe tools that include pandas are:

- [DuckDB (former H2O.ai) benchmarks](https://duckdblabs.github.io/db-benchmark/)
- [TPCH benchmarks](https://pola.rs/posts/benchmarks/)
