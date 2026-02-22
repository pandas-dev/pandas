import re
import subprocess
import sys
import textwrap

from ._base import _get_first_attr
from .time import TimeBenchmark


class _SeparateProcessTimer:
    """
    This class provides a timer that runs a given function in a separate Python
    process.

    The function should return the statement to be timed. This statement is
    executed using the Python timeit module in a new Python process. The
    execution time is then returned.

    #### Attributes
    **subprocess_tmpl** (`str`)
    : The template Python code to be run in the subprocess. It imports necessary
    modules and prints the execution time of the statement.

    **func** (`callable`)
    : The function to be timed. This function should return a string of Python
    code to be executed, or a tuple of two strings: the code to be executed and
    the setup code to be run before timing.

    #### Methods
    **timeit(number)**
    : Run the function's code `number` times in a separate Python process, and
    return the execution time.
    """

    subprocess_tmpl = textwrap.dedent(
        '''
        from __future__ import print_function
        from timeit import timeit, default_timer as timer
        print(repr(timeit(stmt="""{stmt}""", setup="""{setup}""",
                    number={number}, timer=timer)))
    '''
    ).strip()

    def __init__(self, func):
        self.func = func

    def timeit(self, number):
        """
        Run the function's code `number` times in a separate Python process, and
        return the execution time.

        #### Parameters
        **number** (`int`)
        : The number of times to execute the function's code.

        #### Returns
        **time** (`float`)
        : The time it took to execute the function's code `number` times.

        #### Notes
        The function's code is executed in a separate Python process to avoid
        interference from the parent process. The function can return either a
        single string of code to be executed, or a tuple of two strings: the
        code to be executed and the setup code to be run before timing.
        """
        stmt = self.func()
        if isinstance(stmt, tuple):
            stmt, setup = stmt
        else:
            setup = ""
        stmt = textwrap.dedent(stmt)
        setup = textwrap.dedent(setup)
        stmt = stmt.replace(r'"""', r"\"\"\"")
        setup = setup.replace(r'"""', r"\"\"\"")

        code = self.subprocess_tmpl.format(stmt=stmt, setup=setup, number=number)

        res = subprocess.check_output([sys.executable, "-c", code])
        return float(res.strip())


class TimerawBenchmark(TimeBenchmark):
    """
    Represents a benchmark for tracking timing benchmarks run once in
    a separate process.

    This class inherits from `TimeBenchmark` and modifies it to run the
    benchmark function in a separate process. This is useful for isolating the
    benchmark from any potential side effects caused by other Python code
    running in the same process.

    #### Attributes
    **name_regex** (`re.Pattern`)
    : The regular expression used to match the names of functions that should be
    considered as raw timing benchmarks.

    **number** (`int`)
    : The number of times to execute the function's code. By default, the
    function's code is executed once.

    #### Methods
    **_load_vars()**
    : Loads variables for the benchmark from the function's attributes or from
    default values.

    **_get_timer(*param)**
    : Returns a timer that runs the benchmark function in a separate process.

    **do_profile(filename=None)**
    : Raises a ValueError. Raw timing benchmarks cannot be profiled.
    """

    name_regex = re.compile("^(Timeraw[A-Z_].+)|(timeraw_.+)$")

    def _load_vars(self):
        """
        Loads variables for the benchmark from the function's attributes or from
        default values.
        """
        TimeBenchmark._load_vars(self)
        self.number = int(_get_first_attr(self._attr_sources, "number", 1))
        del self.timer

    def _get_timer(self, *param):
        """
        Returns a timer that runs the benchmark function in a separate process.

        #### Parameters
        **param** (`tuple`)
        : The parameters to pass to the benchmark function.

        #### Returns
        **timer** (`_SeparateProcessTimer`)
        : A timer that runs the function in a separate process.
        """
        if param:

            def func():
                self.func(*param)

        else:
            func = self.func
        return _SeparateProcessTimer(func)

    def do_profile(self, filename=None):
        """
        Raises a ValueError. Raw timing benchmarks cannot be profiled.

        #### Parameters
        **filename** (`str`, optional)
        : The name of the file to which to save the profile. Default is None.

        #### Raises
        **ValueError**
        : Always. Raw timing benchmarks cannot be profiled.
        """
        raise ValueError("Raw timing benchmarks cannot be profiled")


export_as_benchmark = [TimerawBenchmark]
