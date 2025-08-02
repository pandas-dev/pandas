import functools
import inspect


class SkipNotImplemented(NotImplementedError):
    """
    Exception raised to indicate a skipped benchmark.

    This exception inherits from `NotImplementedError`. It's used within an ASV
    benchmark to skip the current benchmark for certain parameters or conditions
    that are not implemented or do not apply.

    #### Attributes
    **message** (`str`)
    : A string that provides a more detailed explanation of the skip reason.

    #### Warning
    Use of `SkipNotImplemented` is less efficient than the `@skip_for_params`
    decorator as the setup for the benchmarks and the benchmarks themselves are
    run before the error is raised, thus consuming unnecessary resources. Use
    `@skip_for_params` where possible to avoid running the benchmarks that
    should be skipped.

    #### Notes
    This is mainly provided for backwards compatibility with the behavior of asv
    before 0.5 wherein individual benchmarks could raise and be skipped. From
    0.5 onwards, only the setup function is meant to raise `NotImplemented` for
    skipping parameter sets.

    #### Example
    This exception might be used in a scenario where a benchmark should be
    skipped for certain conditions or parameters:

    ```{code-block} python
    class Simple:
        params = ([False, True])
        param_names = ["ok"]

        def time_failure(self, ok):
            if ok:
                x = 34.2**4.2
            else:
                raise SkipNotImplemented
    ```
    """

    def __init__(self, message=""):
        """
        Initialize a new instance of `SkipNotImplemented`.

        #### Parameters
        **message** (`str`)
        : A string that provides a more detailed explanation of the skip reason.
          Optional; if not provided, defaults to an empty string.
        """
        self.message = message
        super().__init__(self.message)


def skip_for_params(skip_params_list):
    """
    Decorator to set skip parameters for a benchmark function.

    #### Parameters
    **skip_params_list** (`list`):
    A list of tuples, each specifying a combination of parameter values that
    should cause the benchmark function to be skipped.

    #### Returns
    **decorator** (`function`):
    A decorator function that sets the skip parameters for the benchmark
    function.

    #### Notes
    The `skip_for_params` decorator can be used to specify conditions under
    which a benchmark function should be skipped. Each tuple in the list
    represents a combination of parameter values which, if received by the
    benchmark function, will cause that function to be skipped during the
    benchmarking process.

    The decorated function's `skip_params` attribute will be set with the
    provided skip parameters, which will be used during the benchmarking
    process.

    Using this decorator is always more efficient than raising a
    `SkipNotImplemented` exception within the benchmark function, as the
    function setup and execution can be avoided entirely for skipped parameters.

    #### Example
    ```{code-block} python
    class Simple:
        params = ([False, True])
        param_names = ["ok"]

        @skip_for_params([(False, )])
        def time_failure(self, ok):
            if ok:
                x = 34.2**4.2
    ```
    """

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)

        setattr(wrapper, "skip_params", skip_params_list)
        return wrapper

    return decorator


def skip_benchmark(func):
    """
    Decorator to mark a function as skipped for benchmarking.

    #### Parameters
    **func** (function)
    : The function to be marked as skipped.

    #### Returns
    **wrapper** (function)
    : A wrapped function that is marked to be skipped for benchmarking.

    #### Notes
    The `skip_benchmark` decorator can be used to mark a specific function as
    skipped for benchmarking. When the decorated function is encountered during
    benchmarking, it will be skipped and not included in the benchmarking
    process.
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    setattr(wrapper, "skip_benchmark", True)
    return wrapper


def skip_benchmark_if(condition):
    """
    Decorator to skip benchmarking of a function if a condition is met.

    #### Parameters
    **condition** (`bool`)
    : A boolean that indicates whether to skip benchmarking. If `True`,
      the decorated function will be skipped for benchmarking. If `False`,
      the decorated function will be benchmarked as usual.

    #### Returns
    **decorator** (function)
    : A decorator function that sets the condition under which the decorated function
      will be skipped for benchmarking.

    #### Notes
    The `skip_if` decorator can be used to skip the benchmarking of a specific
    function if a condition is met. It is faster than raising
    `SkipNotImplemented` as it skips the `setup()` as well.
    """

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)

        if condition:
            setattr(wrapper, "skip_benchmark", True)
        return wrapper

    return decorator


def skip_params_if(skip_params_list, condition):
    """
    Decorator to set skip parameters for a benchmark function if a condition is met.

    #### Parameters
    **skip_params_list** (`list`):
    A list specifying the skip parameters for the benchmark function.

    **condition** (`bool`)
    : A boolean that indicates whether to set the skip parameters. If `True`,
      the skip parameters will be set for the decorated function. If `False`,
      no parameters will be skipped.

    #### Returns
    **decorator** (function):
    A decorator function that sets the skip parameters for the benchmark function
      if the condition is met.

    #### Notes
    The `skip_params_if` decorator can be used to specify skip parameters for a
      benchmark function if a condition is met.
    """

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)

        if condition:
            setattr(wrapper, "skip_params", skip_params_list)
        return wrapper

    return decorator


def parameterize_class_with(param_dict):
    """
    Class Decorator to set benchmark parameters for a class.

    #### Parameters
    **param_dict** (`dict`):
    A dictionary specifying the parameters for the benchmark class.
    The keys represent the parameter names, and the values are lists
    of values for those parameters.

    #### Returns
    **decorator** (function):
    A class decorator that sets the parameters for the benchmark functions.

    #### Notes
    The `parameterize_class_with` decorator can be used to specify parameters for a
    benchmark class. The parameters are defined as a dictionary, where keys are
    the parameter names and values are lists of respective values. The decorated
    class's `params` and `param_names` attributes will be set with the provided
    parameters and names, which will be used during the benchmarking process.
    This decorator will overwrite any existing `params` and `param_names`
    attributes in the class.
    """

    def decorator(cls):
        if not inspect.isclass(cls):
            raise TypeError(
                "The parameterize_class_with decorator can only be used with classes"
            )
        # Handle the single parameter case separately.
        if len(param_dict) > 1:
            cls.params = list(param_dict.values())
        else:
            cls.params = list(param_dict.values())[0]
        cls.param_names = list(param_dict.keys())
        return cls

    return decorator


def parameterize_func_with(param_dict):
    """
    Function Decorator to set benchmark parameters for a function.

    #### Parameters
    **param_dict** (`dict`):
    A dictionary specifying the parameters for the benchmark function.
    The keys represent the parameter names, and the values are lists
    of values for those parameters.

    #### Returns
    **decorator** (function):
    A function decorator that sets the parameters for the benchmark function.

    #### Notes
    The `parameterize_func_with` decorator can be used to specify parameters for a
    benchmark function. The parameters are defined as a dictionary, where keys are
    the parameter names and values are lists of respective values. The decorated
    function's `params` and `param_names` attributes will be set with the provided
    parameters and names, which will be used during the benchmarking process.
    This decorator will overwrite any existing `params` and `param_names`
    attributes in the function, and it should not be used with methods of a class.
    """

    def decorator(func):
        if inspect.isclass(func) or inspect.ismethod(func):
            raise TypeError(
                "The parameterize_func_with decorator can only be used with functions"
            )
        if len(param_dict) > 1:
            func.params = list(param_dict.values())
        else:
            func.params = list(param_dict.values())[0]
        func.param_names = list(param_dict.keys())
        return func

    return decorator


def parameterize(param_dict):
    """
    Decorator to set benchmark parameters for a function or a class.

    #### Parameters
    **param_dict** (`dict`):
    A dictionary specifying the parameters for the benchmark.
    The keys represent the parameter names, and the values are lists
    of values for those parameters.

    #### Returns
    **decorator** (function):
    A function or class decorator that sets the parameters for the benchmark.

    #### Notes
    The `parameterize` decorator can be used to specify parameters for a
    benchmark function or class. The parameters are defined as a dictionary,
    where keys are the parameter names and values are lists of respective values.
    The decorated function or class's `params` and `param_names` attributes
    will be set with the provided parameters and names, which will be used
    during the benchmarking process.
    """

    def decorator(obj):
        if inspect.isclass(obj):
            return parameterize_class_with(param_dict)(obj)
        elif callable(obj):
            return parameterize_func_with(param_dict)(obj)
        else:
            raise TypeError(
                "The parameterize decorator can only be used with functions or classes"
            )

    return decorator


def timeout_class_at(seconds):
    """
    Class Decorator to set timeout for a class.

    #### Parameters
    **seconds** (`float`)
    : The number of seconds after which the class methods should be timed out.

    #### Returns
    **decorator** (function)
    : A class decorator that sets the timeout for the class.

    #### Notes
    The `timeout_class_at` decorator can be used to specify a timeout for all
    methods in a class. The timeout is stored as an attribute on the class and
    applies to all its methods. Individual methods can override this timeout by
    using the `timeout_func_at` or `timeout_at` decorators.
    """

    def decorator(cls):
        if not inspect.isclass(cls):
            raise TypeError(
                "The timeout_class_with decorator can only be used with classes"
            )
        cls.timeout = seconds
        return cls

    return decorator


def timeout_func_at(seconds):
    """
    Function Decorator to set timeout for a function.

    #### Parameters
    **seconds** (`float`)
    : The number of seconds after which the function should be timed out.

    #### Returns
    **decorator** (function)
    : A function decorator that sets the timeout for the function.

    #### Notes
    The `timeout_func_at` decorator can be used to specify a timeout for a
    specific function. This is particularly useful for benchmarking, where you
    might want to stop execution of functions that take too long. The timeout is
    stored as an attribute on the function.
    """

    def decorator(func):
        if inspect.isclass(func):
            raise TypeError(
                "The timeout_func_with decorator can only be used with functions"
            )

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)

        setattr(wrapper, "timeout", seconds)
        return wrapper

    return decorator


def timeout_at(seconds):
    """
    Decorator to set a timeout for a function or a class.

    #### Parameters
    **seconds** (`float`)
    : The number of seconds after which the function or the class methods should
    be timed out.

    #### Returns
    **decorator** (function)
    : A decorator that sets the timeout for the function or the class.

    #### Notes
    The `timeout_at` decorator can be used to set a specific timeout for a
    function or all methods in a class. If applied to a class, the timeout is
    stored as an attribute on the class and applies to all its methods.
    Individual methods can override this timeout by using the `timeout_func_at`
    or `timeout_at` decorators. If applied to a function, the timeout is stored
    directly on the function.
    """

    def decorator(obj):
        if inspect.isclass(obj):
            return timeout_class_at(seconds)(obj)
        elif callable(obj):
            return timeout_func_at(seconds)(obj)
        else:
            raise TypeError(
                "The parameterize decorator can only be used with functions or classes"
            )

    return decorator


__all__ = [
    "parameterize",
    "skip_benchmark",
    "skip_benchmark_if",
    "skip_for_params",
    "skip_params_if",
    "timeout_at",
]
