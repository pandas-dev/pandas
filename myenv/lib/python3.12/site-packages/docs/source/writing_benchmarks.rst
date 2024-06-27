.. _writing-benchmarks:

Writing benchmarks
==================


.. note::

    The `asv_samples repository
    <https://github.com/airspeed-velocity/asv_samples>`_ has complete examples
    of benchmarks along with continuous integration and can serve as a reference
    for writing and working with benchmarks.

Benchmarks are stored in a Python package, i.e. collection of ``.py``
files in the benchmark suite's ``benchmark`` directory (as defined by
``benchmark_dir`` in the ``asv.conf.json`` file).  The package may
contain arbitrarily nested subpackages, contents of which will also be
used, regardless of the file names.

Within each ``.py`` file, each benchmark is a function or method.  The
name of the function must have a special prefix, depending on the type
of benchmark.  ``asv`` understands how to handle the prefix in either
``CamelCase`` or lowercase with underscores.  For example, to create a
timing benchmark, the following are equivalent::

    def time_range():
        for i in range(1000):
            pass

    def TimeRange():
        for i in range(1000):
            pass

Benchmarks may be organized into methods of classes if desired::

    class Suite:
        def time_range(self):
            for i in range(1000):
                pass

        def time_xrange(self):
            for i in xrange(1000):
                pass

Running benchmarks during development
-------------------------------------

There are some options to ``asv run`` that may be useful when writing
benchmarks.

You may find that ``asv run`` spends a lot of time setting up the
environment each time.  You can have ``asv run`` use an existing
Python environment that already has the benchmarked project and all of
its dependencies installed.  Use the ``--python`` argument to specify
a Python environment to use::

       asv run --python=python

If you don't care about getting accurate timings, but just want to
ensure the code is running, you can add the ``--quick`` argument,
which will run each benchmark only once::

       asv run --quick

In order to display the standard error output (this includes exception tracebacks)
that your benchmarks may produce, pass the ``--show-stderr`` flag::

       asv run --show-stderr

Finally, a quick way to test out the benchmark suite before doing a full run is to use all of
these features together with::

       asv run --python=same --quick --show-stderr --dry-run

.. versionchanged:: 0.6.0
   This replaces the now removed ``asv dev`` command.

You may also want to only do a basic check whether the benchmark suite
is well-formatted, without actually running any benchmarks::

       asv check --python=same

.. _setup-and-teardown:

Setup and teardown functions
----------------------------

If initialization needs to be performed that should not be included in
the timing of the benchmark, include that code in a ``setup`` method
on the class, or add an attribute called ``setup`` to a free function.

For example::

    class Suite:
        def setup(self):
            # load data from a file
            with open("/usr/share/words.txt", "r") as fd:
                self.words = fd.readlines()

        def time_upper(self):
            for word in self.words:
                word.upper()

    # or equivalently...

    words = []
    def my_setup():
        global words
        with open("/usr/share/words.txt", "r") as fd:
            words = fd.readlines()

    def time_upper():
        for word in words:
            word.upper()
    time_upper.setup = my_setup

You can also include a module-level ``setup`` function, which will be
run for every benchmark within the module, prior to any ``setup``
assigned specifically to each function.

Similarly, benchmarks can also have a ``teardown`` function that is
run after the benchmark.  This is useful if, for example, you need to
clean up any changes made to the filesystem.

Note that although different benchmarks run in separate processes, for
a given benchmark repeated measurement (cf. ``repeat`` attribute) and
profiling occur within the same process.  For these cases, the setup
and teardown routines are run multiple times in the same process.

If ``setup`` raises a ``NotImplementedError``, the benchmark is marked
as skipped.

.. note::

   For ``asv`` versions before 0.5 it was possible to raise
   ``NotImplementedError`` from any existing benchmark during its execution and
   the benchmark would be marked as skipped. This behavior was deprecated from
   0.5 onwards.

   .. versionchanged:: 0.6.0

      To keep compatibility with earlier versions, it is possible
      to raise ``asv_runner.benchmark.mark.SkipNotImplemented`` anywhere within a
      Benchmark, though users are advised to use the skip decorators instead as
      they are faster and do not execute the ``setup`` function. See
      :ref:`skipping-benchmarks` for more details.

The ``setup`` method is run multiple times, for each benchmark and for
each repeat.  If the ``setup`` is especially expensive, the
``setup_cache`` method may be used instead, which only performs the
setup calculation once and then caches the result to disk.  It is run
only once also for repeated benchmarks and profiling, unlike
``setup``.  ``setup_cache`` can persist the data for the benchmarks it
applies to in two ways:

- Returning a data structure, which ``asv`` pickles to disk, and
  then loads and passes it as the first argument to each benchmark.

- Saving files to the current working directory (which is a
  temporary directory managed by ``asv``) which are then explicitly
  loaded in each benchmark process.  It is probably best to load
  the data in a ``setup`` method so the loading time is not
  included in the timing of the benchmark.

A separate cache is used for each environment and each commit of the
project being tested and is thrown out between benchmark runs.

For example, caching data in a pickle::

    class Suite:
        def setup_cache(self):
            fib = [1, 1]
            for i in range(100):
                fib.append(fib[-2] + fib[-1])
            return fib

        def track_fib(self, fib):
            return fib[-1]

As another example, explicitly saving data in a file::

    class Suite:
        def setup_cache(self):
            with open("test.dat", "wb") as fd:
                for i in range(100):
                    fd.write('{0}\n'.format(i))

        def setup(self):
            with open("test.dat", "rb") as fd:
                self.data = [int(x) for x in fd.readlines()]

        def track_numbers(self):
            return len(self.data)

The ``setup_cache`` timeout can be specified by setting the
``.timeout`` attribute of the ``setup_cache`` function. The default
value is the maximum of the timeouts of the benchmarks using it.

.. note::

    .. versionchanged:: 0.6.0

        The configuration option ``default_benchmark_timeout``
        can also be set for a project-wide timeout.

.. _benchmark-attributes:

Benchmark attributes
--------------------

Each benchmark can have a number of arbitrary attributes assigned to
it.  The attributes that ``asv`` understands depends on the type of
benchmark and are defined below.  For free functions, just assign the
attribute to the function.  For methods, include the attribute at the
class level.  For example, the following are equivalent::

    def time_range():
        for i in range(1000):
            pass
    time_range.timeout = 120.0

    class Suite:
        timeout = 120.0

        def time_range(self):
            for i in range(1000):
                pass

For the list of attributes, see :doc:`benchmarks`.

.. _parametrized-benchmarks:

Parameterized benchmarks
------------------------

You might want to run a single benchmark for multiple values of some
parameter. This can be done by adding a ``params`` attribute to the
benchmark object::

    def time_range(n):
       for i in range(n):
           pass
    time_range.params = [0, 10, 20, 30]

This will also make the setup and teardown functions parameterized::

    class Suite:
        params = [0, 10, 20]

        def setup(self, n):
            self.obj = range(n)

        def teardown(self, n):
            del self.obj

        def time_range_iter(self, n):
            for i in self.obj:
                pass

If ``setup`` raises a ``NotImplementedError``, the benchmark is marked
as skipped for the parameter values in question.

The parameter values can be any Python objects. However, it is often
best to use only strings or numbers, because these have simple
unambiguous text representations. In the event the ``repr()`` output
is non-unique, the representations will be made unique by suffixing
an integer identifier corresponding to the order of appearance.

When you have multiple parameters, the test is run for all
of their combinations::

     def time_ranges(n, func_name):
         f = {'range': range, 'arange': numpy.arange}[func_name]
         for i in f(n):
             pass

     time_ranges.params = ([10, 1000], ['range', 'arange'])

The test will be run for parameters ``(10, 'range'), (10, 'arange'),
(1000, 'range'), (1000, 'arange')``.

You can also provide informative names for the parameters::

     time_ranges.param_names = ['n', 'function']

These will appear in the test output; if not provided you get default
names such as "param1", "param2".

Note that ``setup_cache`` is not parameterized.

.. _skipping-benchmarks:

Skipping benchmarks
------------------------

.. note::

  This section is only applicable from version ``0.6.0`` on-wards

Conversely, it is possible (typically due to high setup times) that one might
want to skip some benchmarks all-together, or just for some sets of parameters.
This is accomplished by an attribute ``skip_params``, which can be used with the
decorator ``@skip_for_params`` as::

     from asv_runner.benchmarks.mark import skip_for_params
     @skip_for_params([(10, 'arange'), (1000, 'range')])
     def time_ranges(n, func_name):
         f = {'range': range, 'arange': np.arange}[func_name]
         for i in f(n):
             pass

Benchmarks may also be conditionally skipped based on a boolean with ``@skip_benchmark_if``::

     from asv_runner.benchmarks.mark import skip_benchmark_if
     import datetime

     # Skip if not before midday
     @skip_benchmark_if(
         datetime.datetime.now(datetime.timezone.utc).hour >= 12
     )
     def time_ranges(n, func_name):
         f = {'range': range, 'arange': np.arange}[func_name]
         for i in f(n):
             pass

Similarly, for parameters we have ``@skip_params_if``::


     from asv_runner.benchmarks.mark import skip_params_if
     import datetime

     class TimeSuite:
         params = [100, 200, 300, 400, 500]
         param_names = ["size"]

         def setup(self, size):
             self.d = {}
             for x in range(size):
                 self.d[x] = None

         # Skip benchmarking when size is either 100 or 200
         # and the current hour is 12 or later.
         @skip_params_if(
             [(100,), (200,)],
             datetime.datetime.now(datetime.timezone.utc).hour >= 12
         )
         def time_dict_update(self, size):
             d = self.d
             for i in range(size):
                 d[i] = i

.. warning::

   The skips discussed so far, using the decorators will ignore both the
   benchmark, and the ``setup`` function, however, ``setup_cache`` will not be
   affected.

If the onus of preparing the exact parameter sets for ``skip_for_params`` is too
complicated and the ``setup`` function is not too expensive, or if a benchmark
needs to be skipped conditionally but ``skip_*_if`` are not the right choice, there
is also the ``SkipNotImplemented`` exception which can be raised anywhere during
a benchmark run for it to be marked as skipped (``n/a`` in the output table).
This may be used as::

     from asv_runner.benchmarks.mark import SkipNotImplemented
     class SimpleSlow:
         params = ([False, True])
         param_names = ["ok"]
         def time_failure(self, ok):
             if ok:
                 x = 34.2**4.2
             else:
                 raise SkipNotImplemented(f"{ok} is skipped")

Benchmark types
---------------

.. _timing-benchmarks:

Timing
``````

Timing benchmarks have the prefix ``time``.

How ASV runs benchmarks is as follows (pseudocode for main idea)::

     for round in range(`rounds`):
        for benchmark in benchmarks:
            with new process:
                <calibrate `number` if not manually set>
                for j in range(`repeat`):
                    <setup `benchmark`>
                    sample = timing_function(<run benchmark `number` times>) / `number`
                    <teardown `benchmark`>

where the actual ``rounds``, ``repeat``, and ``number`` are :doc:`attributes
of the benchmark <benchmarks>`.

The default timing function is `timeit.default_timer`, which uses the
highest resolution clock available on a given platform to measure the
elapsed wall time. This has the consequence of being more susceptible
to noise from other processes, but the increase in resolution is more
significant for shorter duration tests (particularly on Windows).

Process timing is provided by the function `time.process_time` (POSIX
``CLOCK_PROCESS_CPUTIME``), which measures the CPU time used only by
the current process.  You can change the timer by setting the
benchmark's ``timer`` attribute, for example to `time.process_time`
to measure process time.

.. note::

   One consequence of using `time.process_time` is that the time
   spent in child processes of the benchmark is not included.
   Multithreaded benchmarks also return the total CPU time
   counting all CPUs. In these cases you may want to measure the
   wall clock time, by setting the
   ``timer = timeit.default_timer`` benchmark attribute.

For best results, the benchmark function should contain as little as
possible, with as much extraneous setup moved to a ``setup`` function::

    class Suite:
        def setup(self):
            # load data from a file
            with open("/usr/share/words.txt", "r") as fd:
                self.words = fd.readlines()

        def time_upper(self):
            for word in self.words:
                word.upper()

How ``setup`` and ``teardown`` behave for timing benchmarks
is similar to the Python ``timeit`` module, and the behavior is controlled
by the ``number`` and ``repeat`` attributes.

For the list of benchmark attributes, see :doc:`benchmarks`.

.. _memory-benchmarks:

Memory
``````

Memory benchmarks have the prefix ``mem``.

Memory benchmarks track the size of Python objects.  To write a memory
benchmark, write a function that returns the object you want to track::

    def mem_list():
        return [0] * 256

The `asizeof <http://pythonhosted.org/Pympler/asizeof.html>`__ module
is used to determine the size of Python objects.  Since ``asizeof``
includes the memory of all of an object's dependencies (including the
modules in which their classes are defined), a memory benchmark
instead calculates the incremental memory of a copy of the object,
which in most cases is probably a more useful indicator of how much
space *each additional* object will use.  If you need to do something
more specific, a generic :ref:`tracking` benchmark can be used
instead.

For details, see :doc:`benchmarks`.

.. note::

    The memory benchmarking feature is still experimental.
    ``asizeof`` may not be the most appropriate metric to use.

.. note::

    The memory benchmarks are not supported on PyPy.

.. _peak-memory:

Peak Memory
```````````

Peak memory benchmarks have the prefix ``peakmem``.

Peak memory benchmark tracks the maximum resident size (in bytes) of
the process in memory. This does not necessarily count memory paged
on-disk, or that used by memory-mapped files.  To write a peak memory
benchmark, write a function that does the operation whose maximum
memory usage you want to track::

    def peakmem_list():
        [0] * 165536


.. note::

   The peak memory benchmark also counts memory usage during the
   ``setup`` routine, which may confound the benchmark results. One
   way to avoid this is to use ``setup_cache`` instead.

For details, see :doc:`benchmarks`.


.. _raw-timing-benchmarks:

Raw timing benchmarks
`````````````````````

For some timing benchmarks, for example measuring the time it takes to
import a module, it is important that they are run separately in a new
Python process.

Measuring execution time for benchmarks run once in a new Python process
can be done with ``timeraw_*`` timing benchmarks::

    def timeraw_import_inspect():
        return """
        import inspect
        """

Note that these benchmark functions should return a string,
corresponding to the code that will be run.

Importing a module takes a meaningful amount of time only the first time
it is executed, therefore a fresh interpreter is used for each iteration of
the benchmark. The string returned by the benchmark function is executed in a
subprocess.

Note that the setup and setup_cache are performed in the base benchmark
process, so that the setup done by them is not available in the benchmark code.
To perform setup also in the benchmark itself, you can return a second string:

    def timeraw_import_inspect():
        code = "import inspect"
        setup = "import ast"
        return code, setup

The raw timing benchmarks have the same parameters as ordinary timing benchmarks,
but ``number`` is by default 1, and ``timer`` is ignored.

.. note::

   Timing standard library modules is possible as long as they are not
   `built-in`_ or brought in by importing the ``timeit`` module (which
   further imports ``gc``, ``sys``, ``time``, and ``itertools``).

.. _built-in: https://hg.python.org/cpython/file/tip/Modules/Setup.dist


Imports
```````

You can use raw timing benchmarks to measure import times.


.. _tracking:

Tracking (Generic)
``````````````````

It is also possible to use ``asv`` to track any arbitrary numerical
value.  "Tracking" benchmarks can be used for this purpose and use the
prefix ``track``.  These functions simply need to return a numeric
value.  For example, to track the number of objects known to the
garbage collector at a given state::

    import gc

    def track_num_objects():
        return len(gc.get_objects())
    track_num_objects.unit = "objects"

For details, see :doc:`benchmarks`.


Benchmark versioning
--------------------

When you edit benchmark's code in the benchmark suite, this often
changes what is measured, and previously measured results should be
discarded.

Airspeed Velocity records with each benchmark measurement a "version
number" for the benchmark. By default, it is computed by hashing the
benchmark source code text, including any ``setup`` and
``setup_cache`` routines.  If there are changes in the source code of
the benchmark in the benchmark suite, the version number changes, and
``asv`` will ignore results whose version number is different from the
current one.

It is also possible to control the versioning of benchmark results
manually, by setting the ``.version`` attribute for the benchmark. The
version number, i.e. content of the attribute, can be any Python
string. ``asv`` only checks whether the version recorded with a
measurement matches the current version, so you can use any versioning
scheme.

See :doc:`benchmarks` for reference documentation.
