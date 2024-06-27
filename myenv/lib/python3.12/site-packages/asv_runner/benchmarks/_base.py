import cProfile as profile
import inspect
import itertools
import math
import os
import re
import textwrap
from collections import Counter
from hashlib import sha256


def _get_attr(source, name, ignore_case=False):
    """
    Retrieves an attribute from a source by its name.

    #### Parameters
    **source** (`object`)
    : The source from which to get the attribute.

    **name** (`str`)
    : The name of the attribute.

    **ignore_case** (`bool`, optional)
    : Whether to ignore case when comparing attribute names. Defaults to `False`.

    #### Returns
    **attr** (`object` or `None`)
    : The attribute if it is found, else `None`.

    #### Raises
    **ValueError**
    : If more than one attribute with the given name exists and `ignore_case` is `True`.
    """
    if not ignore_case:
        return getattr(source, name, None)
    attrs = [getattr(source, key) for key in dir(source) if key.lower() == name.lower()]

    if len(attrs) > 1:
        raise ValueError(f"{source.__name__} contains multiple {name} functions.")
    elif len(attrs) == 1:
        return attrs[0]
    else:
        return None


def _get_all_attrs(sources, name, ignore_case=False):
    """
    Yields attributes from a list of sources by their name.

    #### Parameters
    **sources** (`List[object]`)
    : The list of sources from which to get the attribute.

    **name** (`str`)
    : The name of the attribute.

    **ignore_case** (`bool`, optional)
    : Whether to ignore case when comparing attribute names. Defaults to `False`.

    #### Yields
    **val** (`object`)
    : The attribute if it is found in the source.
    """
    for source in sources:
        val = _get_attr(source, name, ignore_case=ignore_case)
        if val is not None:
            yield val


def _get_first_attr(sources, name, default, ignore_case=False):
    """
    Retrieves the first attribute from a list of sources by its name.

    #### Parameters
    **sources** (`List[object]`)
    : The list of sources from which to get the attribute.

    **name** (`str`)
    : The name of the attribute.

    **default** (`object`)
    : The default value to return if no attribute is found.

    **ignore_case** (`bool`, optional)
    : Whether to ignore case when comparing attribute names. Defaults to `False`.

    #### Returns
    **attr** (`object`)
    : The first attribute found or the default value if no attribute is found.
    """
    for val in _get_all_attrs(sources, name, ignore_case=ignore_case):
        return val
    return default


def get_setup_cache_key(func):
    """
    Retrieves the cache key for a function's setup.

    #### Parameters
    **func** (`function`)
    : The function for which to get the cache key.

    #### Returns
    **cache_key** (`str` or `None`)
    : The cache key if the function is not `None`, else `None`.

    #### Notes
    The cache key is a string composed of the function's module name and the line
    number where the function's source code starts.
    """
    if func is None:
        return None

    module = inspect.getmodule(func)
    mname = ".".join(module.__name__.split(".", 1)[1:])
    if not mname:
        mname = inspect.getsourcefile(func)

    return f"{mname}:{inspect.getsourcelines(func)[1]}"


def get_source_code(items):
    """
    Extracts, concatenates, and dedents the source code of the given items.

    #### Parameters
    **items** (`Iterable[object]`)
    : An iterable of items, typically functions or methods, for which to extract the
    source code.

    #### Returns
    **source_code** (`str`)
    : The concatenated and dedented source code of the items.

    #### Notes
    The function retrieves the source code of each item. If the item has a
    `pretty_source` attribute, it uses that as the source code. Otherwise, it
    attempts to use the `inspect` module's `getsourcelines` function to extract
    the source code.

    The function also adds class names to methods and properly indents the
    source code.  If the source code belongs to a method, the function retrieves
    the class name and prepends it to the source code, properly indenting it to
    reflect its position within the class. If the source code belongs to the
    same class as the previous item, only the indentation is adjusted.
    """
    sources = []
    prev_class_name = None

    for func in items:
        # custom source
        if hasattr(func, "pretty_source"):
            src = textwrap.dedent(func.pretty_source).lstrip()
        # original source
        else:
            try:
                lines, _ = inspect.getsourcelines(func)
            except TypeError:
                continue

            if not lines:
                continue

            src = "\n".join(line.rstrip() for line in lines)
            src = textwrap.dedent(src)

        class_name = None
        if inspect.ismethod(func):
            # Add class name
            if hasattr(func, "im_class"):
                class_name = func.im_class.__name__
            elif hasattr(func, "__qualname__"):
                names = func.__qualname__.split(".")
                if len(names) > 1:
                    class_name = names[-2]

        if class_name and prev_class_name != class_name:
            src = "class {}:\n    {}".format(class_name, src.replace("\n", "\n    "))
        elif class_name:
            src = "    " + src.replace("\n", "\n    ")

        sources.append(src)
        prev_class_name = class_name

    return "\n\n".join(sources).rstrip()


def _get_sourceline_info(obj, basedir):
    """
    Retrieves the source file and line number information of the given object.

    #### Parameters
    **obj** (`object`)
    : The object for which to retrieve source file and line number information. This is
    typically a function or a method.

    **basedir** (`str`)
    : The base directory relative to which the source file path should be expressed.

    #### Returns
    **sourceline_info** (`str`)
    : A string containing the relative path of the source file and the line number where
    the object is defined, in the format `' in {filename}:{lineno}'`. If the source file
    or line number cannot be determined, an empty string is returned.

    #### Notes
    The function uses the `inspect` module's `getsourcefile` and
    `getsourcelines` functions to determine the source file and line number of
    the object, respectively.  The source file path is converted to a path
    relative to `basedir` using `os.path.relpath`.
    """
    try:
        fn = inspect.getsourcefile(obj)
        fn = os.path.relpath(fn, basedir)
        _, lineno = inspect.getsourcelines(obj)
        return f" in {fn !s}:{lineno !s}"
    except Exception:
        return ""


def _check_num_args(root, benchmark_name, func, min_num_args, max_num_args=None):
    """
    Verifies if the function under benchmarking accepts a correct number of arguments.

    #### Parameters
    **root** (`str`)
    : The root directory for the function's source file (used to print detailed error
    messages).

    **benchmark_name** (`str`)
    : The name of the benchmark for which the function is being checked (used in error
    messages).

    **func** (`function`)
    : The function to check for correct number of arguments.

    **min_num_args** (`int`)
    : The minimum number of arguments the function should accept.

    **max_num_args** (`int`, optional)
    : The maximum number of arguments the function should accept. If not provided,
    `max_num_args` is assumed to be the same as `min_num_args`.

    #### Returns
    **validity** (`bool`)
    : True if the function accepts a correct number of arguments, False otherwise.

    #### Notes
    The function uses the `inspect` module's `getfullargspec` function to determine the
    number of arguments the function accepts. It correctly handles functions, methods,
    variable argument lists, and functions with default argument values. In case of any
    error or if the function does not accept a correct number of arguments, an error
    message is printed to standard output.
    """
    if max_num_args is None:
        max_num_args = min_num_args

    try:
        info = inspect.getfullargspec(func)
    except Exception as exc:
        print(
            f"{benchmark_name !s}: failed to check "
            f"({func !r}{_get_sourceline_info(func, root) !s}): {exc !s}"
        )
        return True

    max_args = len(info.args)

    if inspect.ismethod(func):
        max_args -= 1

    min_args = max_args if info.defaults is None else max_args - len(info.defaults)
    if info.varargs is not None:
        max_args = math.inf

    ok = (min_args <= max_num_args) and (min_num_args <= max_args)
    if not ok:
        args_str = min_args if min_args == max_args else f"{min_args}-{max_args}"
        if min_num_args == max_num_args:
            num_args_str = min_num_args
        else:
            num_args_str = f"{min_num_args}-{max_num_args}"
        print(
            f"{benchmark_name !s}: wrong number of arguments "
            f"(for {func !r}{_get_sourceline_info(func, root) !s}):",
            f"expected {num_args_str}, " f"has {args_str}",
        )

    return ok


def _repr_no_address(obj):
    """
    Returns a string representing the object, but without its memory address.

    #### Parameters
    **obj** (`object`)
    : The object to represent.

    #### Returns
    **representation** (`str`)
    : A string representation of the object without its memory address.

    #### Notes
    When Python's built-in `repr` function is used on an object, it often includes
    the memory address of the object. In some cases, this might not be desirable
    (for example, when comparing object representations in unit tests, where the
    memory address is not relevant). This function provides a way to get a string
    representation of an object without its memory address.

    The function works by first getting the `repr` of the object, then using a
    regular expression to detect and remove the memory address if it's present.
    To avoid false positives, the function also gets the `repr` of the object
    using the `object` class's `__repr__` method (which always includes the
    address), and only removes the address from the original `repr` if it matches
    the address in the `object.__repr__`.

    Please note, this function is not guaranteed to remove the memory address for
    all objects. It is primarily intended to work for objects that have a `repr`
    similar to the default one provided by the `object` class.
    """
    result = repr(obj)
    address_regex = re.compile(r"^(<.*) at (0x[\da-fA-F]*)(>)$")
    match = address_regex.match(result)
    if match:
        suspected_address = match[2]
        # Double check this is the actual address
        default_result = object.__repr__(obj)
        match2 = address_regex.match(default_result)
        if match2:
            known_address = match2[2]
            if known_address == suspected_address:
                result = match[1] + match[3]

    return result


def _validate_params(params, param_names, name):
    """
    Validates the params and param_names attributes and returns validated lists.

    #### Parameters
    **params** (`list`)
    : List of parameters for the function to be benchmarked.

    **param_names** (`list`)
    : List of names for the parameters.

    **name** (`str`)
    : The name of the benchmark.

    #### Returns
    **params**, **param_names** (`list`, `list`)
    : The validated parameter and parameter name lists.
    """

    try:
        param_names = [str(x) for x in list(param_names)]
    except ValueError:
        raise ValueError(f"{name}.param_names is not a list of strings")

    try:
        params = list(params)
    except ValueError:
        raise ValueError(f"{name}.params is not a list")

    if params and not isinstance(params[0], (tuple, list)):
        params = [params]
    else:
        params = [list(entry) for entry in params]

    if len(param_names) != len(params):
        param_names = param_names[: len(params)]
        param_names += [
            "param%d" % (k + 1,) for k in range(len(param_names), len(params))
        ]

    return params, param_names


def _unique_param_ids(params):
    """
    Processes the params list to handle duplicate names within parameter sets,
    ensuring unique IDs.

    #### Parameters
    **params** (`list`)
    : List of parameters. Each entry is a list representing a set of parameters.

    #### Returns
    **params** (`list`)
    : List of parameters with duplicate names within each set handled.
    If there are duplicate names, they are renamed with a numerical suffix to
    ensure unique IDs.
    """

    params = [[_repr_no_address(item) for item in entry] for entry in params]
    for i, param in enumerate(params):
        if len(param) != len(set(param)):
            counter = Counter(param)
            dupe_dict = {name: 0 for name, count in counter.items() if count > 1}
            for j in range(len(param)):
                name = param[j]
                if name in dupe_dict:
                    param[j] = f"{name} ({dupe_dict[name]})"
                    dupe_dict[name] += 1
            params[i] = param
    return params


class Benchmark:
    """
    Class representing a single benchmark. The class encapsulates
    functions and methods that can be marked as benchmarks, along with
    setup and teardown methods, timing and other configuration.

    #### Notes
    The class uses regex to match method names that will be considered
    as benchmarks. The matched functions are then processed for
    benchmarking using various helper methods.

    By default, a benchmark's timeout is set to 60 seconds.
    """

    # The regex of the name of function or method to be considered as
    # this type of benchmark.  The default in the base class, will
    # match nothing.
    name_regex = re.compile("^$")

    def __init__(self, name, func, attr_sources):
        """
        Initialize a new instance of `Benchmark`.

        #### Parameters
        **name** (`str`)
        : The name of the benchmark.

        **func** (`function`)
        : The function to benchmark.

        **attr_sources** (`list`)
        : List of sources from which attributes of the benchmark will be drawn.
        These attributes include setup, teardown, timeout, etc.

        #### Attributes
        **pretty_name** (`str`)
        : A user-friendly name for the function being benchmarked, if available.

        **_setups** (`list`)
        : List of setup methods to be executed before the benchmark.

        **_teardowns** (`list`)
        : List of teardown methods to be executed after the benchmark.

        **_setup_cache** (`function`)
        : A special setup method that is only run once per parameter set.

        **setup_cache_key** (`str`)
        : A unique key for the setup cache.

        **setup_cache_timeout** (`float`)
        : The time after which the setup cache should be invalidated.

        **timeout** (`float`)
        : The maximum time the benchmark is allowed to run before it is aborted.

        **code** (`str`)
        : The source code of the function to be benchmarked and its setup methods.

        **version** (`str`)
        : A version string derived from a hash of the code.

        **_params** (`list`)
        : List of parameters for the function to be benchmarked.

        **param_names** (`list`)
        : List of names for the parameters.

        **_current_params** (`tuple`)
        : The current set of parameters to be passed to the function during the
        benchmark.

        **params** (`list`)
        : The list of parameters with unique representations for exporting.

        **_skip_tuples** (`list`)
        : List of tuples representing parameter combinations to be skipped
        before calling the setup method.

        #### Raises
        **ValueError**
        : If `param_names` or `_params` is not a list or if the number of
        parameters does not match the number of parameter names.
        """
        self.name = name
        self.func = func
        self.pretty_name = getattr(func, "pretty_name", None)
        self._attr_sources = attr_sources
        self._setups = list(_get_all_attrs(attr_sources, "setup", True))[::-1]
        self._teardowns = list(_get_all_attrs(attr_sources, "teardown", True))
        self._setup_cache = _get_first_attr(attr_sources, "setup_cache", None)
        self.setup_cache_key = get_setup_cache_key(self._setup_cache)
        self.setup_cache_timeout = _get_first_attr([self._setup_cache], "timeout", None)
        self.timeout = _get_first_attr(attr_sources, "timeout", None)
        self.code = get_source_code([self.func] + self._setups + [self._setup_cache])
        code_text = self.code.encode("utf-8")
        code_hash = sha256(code_text).hexdigest()
        self.version = str(_get_first_attr(attr_sources, "version", code_hash))
        self.type = "base"
        self.unit = "unit"

        self._redo_setup_next = False

        self._params = _get_first_attr(attr_sources, "params", [])
        self.param_names = _get_first_attr(attr_sources, "param_names", [])
        self._current_params = ()

        self._params, self.param_names = _validate_params(
            self._params, self.param_names, self.name
        )

        # Fetch skip parameters
        self._skip_tuples = _get_first_attr(attr_sources, "skip_params", [])

        # Exported parameter representations
        self.params = _unique_param_ids(self._params)

    def __repr__(self):
        return f"<{self.__class__.__name__} {self.name}>"

    def set_param_idx(self, param_idx):
        """
        Set the current parameter values for the benchmark based on a parameter
        index.

        This method updates the `_current_params` attribute with the set of
        parameter values that correspond to the provided parameter index.

        #### Parameters
        **param_idx** (`int`)
        : The index of the desired parameter set in the Cartesian product of
        `_params` attribute list.

        #### Raises
        **ValueError**
        : If the provided parameter index is not valid. This could occur if the
        index does not correspond to any element in the Cartesian product of the
        `_params` list.
        """
        try:
            (self._current_params,) = itertools.islice(
                itertools.product(*self._params), param_idx, param_idx + 1
            )
        except ValueError:
            raise ValueError(
                f"Invalid benchmark parameter permutation index: {param_idx!r}"
            )

    def insert_param(self, param):
        """
        Inserts a parameter at the beginning of the current parameter list.

        This method modifies the `_current_params` attribute, inserting the provided
        parameter value at the front of the parameter tuple.

        #### Parameters
        **param** (`Any`)
        : The parameter value to insert at the front of `_current_params`.
        """
        self._current_params = tuple([param] + list(self._current_params))

    def check(self, root):
        """
        Checks call syntax (argument count) for benchmark's setup, call, and teardown.

        #### Parameters
        **root** (`Any`)
        : The root context for checking argument count in setup, call and teardown.

        #### Returns
        **result** (`bool`)
        : `True` if correct argument count is used in all methods, `False` otherwise.

        #### Notes
        The call syntax is checked only based on the number of arguments.  It
        also sets the current parameters for the benchmark if they exist.  The
        number of arguments required by setup, call, and teardown methods may
        increase if a setup cache is defined.
        """
        # Check call syntax (number of arguments only...)
        ok = True

        if self._params:
            self.set_param_idx(0)

        min_num_args = len(self._current_params)
        max_num_args = min_num_args

        if self.setup_cache_key is not None:
            ok = ok and _check_num_args(
                root, f"{self.name}: setup_cache", self._setup_cache, 0
            )
            max_num_args += 1

        for setup in self._setups:
            ok = ok and _check_num_args(
                root, f"{self.name}: setup", setup, min_num_args, max_num_args
            )

        ok = ok and _check_num_args(
            root, f"{self.name}: call", self.func, min_num_args, max_num_args
        )

        for teardown in self._teardowns:
            ok = ok and _check_num_args(
                root,
                f"{self.name}: teardown",
                teardown,
                min_num_args,
                max_num_args,
            )

        return ok

    def do_setup(self):
        if tuple(self._current_params) in self._skip_tuples:
            # Skip
            return True
        try:
            for setup in self._setups:
                setup(*self._current_params)
        except NotImplementedError as e:
            # allow skipping test
            print(f"asv: skipped: {e !r} ")
            return True
        return False

    def redo_setup(self):
        if not self._redo_setup_next:
            self._redo_setup_next = True
            return
        self.do_teardown()
        self.do_setup()

    def do_teardown(self):
        if tuple(self._current_params) in self._skip_tuples:
            # Skip
            return
        for teardown in self._teardowns:
            teardown(*self._current_params)

    def do_setup_cache(self):
        if self._setup_cache is not None:
            return self._setup_cache()

    def do_run(self):
        if tuple(self._current_params) in self._skip_tuples:
            # Skip
            return
        return self.run(*self._current_params)

    def do_profile(self, filename=None):
        """
        Executes the benchmark's function with profiling using `cProfile`.

        #### Parameters
        **filename** (`str`, optional)
        : The name of the file where the profiling data should be saved. If not
          provided, the profiling data will not be saved.

        #### Raises
        **RuntimeError**
        : If the `cProfile` module couldn't be imported.

        #### Notes
        The method uses an inner function `method_caller` to call the function
        to be profiled. The function and its parameters should be available in
        the scope where `method_caller` is called.

        The `cProfile` module should be available, or else a `RuntimeError` is
        raised. If a `filename` is provided, the profiling results will be saved
        to that file.
        """
        if tuple(self._current_params) in self._skip_tuples:
            # Skip
            return

        def method_caller():
            run(*params)  # noqa:F821 undefined name

        if profile is None:
            raise RuntimeError("cProfile could not be imported")

        if filename is not None:
            if hasattr(method_caller, "func_code"):
                code = method_caller.func_code
            else:
                code = method_caller.__code__

            self.redo_setup()

            profile.runctx(
                code, {"run": self.func, "params": self._current_params}, {}, filename
            )
