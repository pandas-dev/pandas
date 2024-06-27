Tuning timing measurements
==========================

The results from timing benchmarks are generally variable.

Performance variations occur on different time scales. For timing
benchmarks repeated immediately after each other, there is always some
jitter in the results, due to operating system scheduling and other
sources.  For timing benchmarks run at more widely separated times,
systematic differences changing on long time scales can appear, for
example from changes in the background system load or built-in CPU
mechanisms for power and heat management.

Airspeed Velocity has mechanisms to deal with these variations.  For
dealing with short-time variations, you can use the ``sample_time``,
``number`` and ``repeat`` attributes of timing benchmarks to control
how results are sampled and averaged.  For long-time variations, you
can use the ``rounds`` attribute and ``--interleave-rounds``,
``--append-samples``, and ``-a rounds=4`` command line options to
run timing benchmarks at more widely spaced times, in order to average
over long-time performance variations.

If you are planning to capture historical benchmark data for most
commits, very accurate timings are not necessary.  The detection of
regressions in historical benchmark data used in ``asv`` is designed
to be statistically robust and tolerates fair amounts of noise.
However, if you are planning to use ``asv continuous`` and ``asv
compare``, accurate results are more important.

Library settings
----------------

If your code uses 3rd party libraries, you may want to check their
settings before benchmarking.  In particular, such libraries may use
automatic multithreading, which may affect runtime performance in
surprising ways.  If you are using libraries such as OpenBLAS, Intel
MKL, or OpenMP, benchmark results may become easier to understand by
forcing single-threaded operation. For these three, this can be
typically done by setting environment variables::

    OPENBLAS_NUM_THREADS=1
    MKL_NUM_THREADS=1
    OMP_NUM_THREADS=1


Tuning machines for benchmarking
--------------------------------

Especially if you are using a laptop computer for which the heat and
power management is an issue, getting reliable results may require too
long averaging times to be practical.

To improve the situation it is possible to optimize the usage and
settings of your machine to minimize the variability in timing
benchmarks.  Generally, while running benchmarks there should not be
other applications actively using CPU, or you can run ``asv`` pinned
to a CPU core not used by other processes.  You should also force the
CPU frequency or power level settings to a fixed value.

The `pyperf <https://pyperf.readthedocs.io/>`__ project has `documentation
on how to tune machines for benchmarking
<https://pyperf.readthedocs.io/en/latest/system.html>`__.  The simplest
way to apply basic tuning on Linux using ``pyperf`` is to run::

    sudo python -mpyperf system tune

This will modify system settings that can be only changed as root, and
you should read the ``pyperf`` documentation on what it precisely does.
This system tuning also improves results for ``asv``.  To achieve CPU
affinity pinning with ``asv`` (e.g. to an isolated CPU), you should
use :ref:`the --cpu-affinity option <cmd-asv-run>`.

It is also useful to note that configuration changes and operating
system upgrades on the benchmarking machine can change the baseline
performance of the machine. For absolutely best results, you may then
want to use a dedicated benchmarking machine that is not used for
anything else. You may also want to carefully select a long-term
supported operating system, such that you can only choose to install
security upgrades.
