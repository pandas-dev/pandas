Benchmark types and attributes
==============================

.. only:: not man

   .. contents::

.. warning::

   .. versionchanged:: 0.6.0

      The code for these have now been moved to be in ``asv_runner``, and the
      rest of the documentation may be outdated.

Benchmark types
---------------

The following benchmark types are recognized:

- ``def time_*()``: measure time taken by the function. See :ref:`timing-benchmarks`.
- ``def timeraw_*()``: measure time taken by the function, after interpreter start. See :ref:`raw-timing-benchmarks`.
- ``def mem_*()``: measure memory size of the object returned.  See :ref:`memory-benchmarks`.
- ``def peakmem_*()``: measure peak memory size of the process when calling the function.
  See :ref:`peak-memory`.
- ``def track_*()``: use the returned numerical value as the benchmark result
  See :ref:`tracking`.

.. note::

   .. versionadded:: 0.6.2

        External benchmarks may be defined through ``asv_runner`` and a list of
        benchmark plugins (like ``asv_bench_memray``) may be found here, with
        samples at `asv_samples
        <https://github.com/airspeed-velocity/asv_samples>`_.

Benchmark attributes
--------------------

Benchmark attributes can either be applied directly to the benchmark function::

    def time_something():
        pass

    time_something.timeout = 123

or appear as class attributes::

    class SomeBenchmarks:
        timeout = 123

        def time_something(self):
            pass

Different benchmark types have their own sets of applicable
attributes.  Moreover, the following attributes are applicable to all
benchmark types:

- ``timeout``: The amount of time, in seconds, to give the benchmark
  to run before forcibly killing it.  Defaults to 60 seconds.

- ``benchmark_name``: If given, used as benchmark function name instead of generated one
  ``<module>.<class>.<function>``.

- ``pretty_name``: If given, used to display the benchmark name instead of the
  benchmark function name.

- ``pretty_source``: If given, used to display a custom version of the benchmark source.

- ``version``: Used to determine when to invalidate old benchmark
  results.  Benchmark results produced with a different value of the
  version than the current value will be ignored.  The value can be
  any Python string (or other object, ``str()`` will be taken).

  Default (if ``version=None`` or not given): hash of the source code
  of the benchmark function and setup and setup_cache methods. If the
  source code of any of these changes, old results become invalidated.

- ``setup``: function to be called as a setup function for the benchmark
  See :ref:`setup-and-teardown` for discussion.

- ``teardown``: function to be called as a teardown function for the benchmark
  See :ref:`setup-and-teardown` for discussion.

- ``setup_cache``: function to be called as a cache setup function.
  See :ref:`setup-and-teardown` for discussion.

- ``param_names``: list of parameter names
  See :ref:`parametrized-benchmarks` for discussion.

- ``params``: list of lists of parameter values.
  If there is only a single parameter, may also be a list of parameter values.
  See :ref:`parametrized-benchmarks` for discussion.

  Example::

     def setup_func(n, func):
         print(n, func)

     def teardown_func(n, func):
         print(n, func)

     def time_ranges(n, func):
         for i in func(n):
             pass

     time_ranges.setup = setup_func
     time_ranges.param_names = ['n', 'func']
     time_ranges.params = ([10, 1000], [range, numpy.arange])

  The benchmark will be run for parameters ``(10, range), (10,
  numpy.arange), (1000, range), (1000, numpy.arange)``. The setup and
  teardown functions will also obtain these parameters.

  Note that ``setup_cache`` is not parameterized.

  For the purposes of identifying benchmarks in the UI, ``repr()`` is called
  on the elements of ``params``. In the event these strings contain memory
  addresses, those adresses are stripped to allow comparison across runs.
  Additionally, if this results in a non-unique mapping, each duplicated
  element will be suffixed with a distinct integer identifier corresponding
  to order of appearance.

Timing benchmarks
`````````````````

- ``warmup_time``: ``asv`` will spend this time (in seconds) in calling
  the benchmarked function repeatedly, before starting to run the actual
  benchmark. If not specified, ``warmup_time`` defaults to 0.1 seconds
  (on PyPy, the default is 1.0 sec).

- ``rounds``: How many rounds to run the benchmark in (default: 2).
  The rounds run different timing benchmarks in an interleaved order,
  allowing to sample over longer periods of background performance
  variations (e.g. CPU power levels).

- ``repeat``: The number measurement samples to collect per round.
  Each sample consists of running the benchmark ``number`` times.
  The median time from all samples collected in all roudns is used
  as the final measurement result.

  ``repeat`` can be a tuple ``(min_repeat, max_repeat, max_time)``.
  In this case, the measurement first collects at least ``min_repeat``
  samples, and continues until either ``max_repeat`` samples are collected
  or the collection time exceeds ``max_time``.

  When not provided (``repeat`` set to 0), the default value is
  ``(1, 10, 20.0)`` if ``rounds==1`` and ``(1, 5, 10.0)`` otherwise.

- ``number``: Manually choose the number of iterations in each sample.
  If ``number`` is specified, ``sample_time`` is ignored.
  Note that ``setup`` and ``teardown`` are not run between iterations:
  ``setup`` runs first, then the timed benchmark routine is called
  ``number`` times, and after that ``teardown`` runs.

- ``sample_time``: ``asv`` will automatically select ``number`` so that
  each sample takes approximatively ``sample_time`` seconds.  If not
  specified, ``sample_time`` defaults to 10 milliseconds.

- ``min_run_count``: the function is run at least this many times during
  benchmark. Default: 2

- ``timer``: The timing function to use, which can be any source of
  monotonically increasing numbers, such as ``time.clock``, ``time.time``
  or ``time.process_time``.  If it's not provided, it defaults to
  ``timeit.default_timer``, but other useful values are
  ``process_time``, for which ``asv`` provides a backported version for
  versions of Python prior to 3.3.

  .. versionchanged:: 0.4

     Previously, the default timer measured process time, which was chosen
     to minimize noise from other processes. However, on Windows, this is
     only available at a resolution of 15.6ms, which is greater than the
     recommended benchmark runtime of 10ms. Therefore, we default to the
     highest resolution clock on any platform.

The ``sample_time``, ``number``, ``repeat``, and ``timer`` attributes
can be adjusted in the ``setup()`` routine, which can be useful for
parameterized benchmarks.


Tracking benchmarks
```````````````````

- ``unit``: The unit of the values returned by the benchmark.  Used
  for display in the web interface.


Environment variables
---------------------

When ``asv`` runs benchmarks, several environment variables are
defined, see :doc:`env_vars`.
