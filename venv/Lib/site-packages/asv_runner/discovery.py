import importlib
import inspect
import json
import os
import pkgutil
import traceback

from ._aux import update_sys_path
from .benchmarks import benchmark_types


def _get_benchmark(attr_name, module, klass, func):
    """
    Retrieves benchmark function based on attribute name, module, class, and
    function.

    #### Parameters
    **attr_name** (`str`)
    : The attribute name of the function.

    **module** (module)
    : The module where the function resides.

    **klass** (class or None)
    : The class defining the function, or None if not applicable.

    **func** (function)
    : The function to be benchmarked.

    #### Returns
    **benchmark** (Benchmark instance or None)
    : A benchmark instance with the name of the benchmark, the function to be
    benchmarked, and its sources. Returns None if no matching benchmark is found
    or the function is marked to be skipped.

    #### Notes
    The function tries to get the `benchmark_name` from `func`. If it fails, it
    uses `attr_name` to match with the name regex in the benchmark types.  If a
    match is found, it creates a new benchmark instance and returns it.  If no
    match is found or the function is marked to be skipped, it returns None.
    """
    # Check if the function has been marked to be skipped
    if getattr(func, "skip_benchmark", False):
        return

    try:
        name = func.benchmark_name
    except AttributeError:
        name = None
        search = attr_name
    else:
        search = name.split(".")[-1]

    for cls in benchmark_types:
        if cls.name_regex.match(search):
            break
    else:
        return
    # relative to benchmark_dir
    mname_parts = module.__name__.split(".", 1)[1:]
    if klass is None:
        if name is None:
            name = ".".join(mname_parts + [func.__name__])
        sources = [func, module]
    else:
        instance = klass()
        func = getattr(instance, attr_name)
        if name is None:
            name = ".".join(mname_parts + [klass.__name__, attr_name])
        sources = [func, instance, module]
    return cls(name, func, sources)


def disc_modules(module_name, ignore_import_errors=False):
    """
    Recursively imports a module and all sub-modules in the package.

    #### Parameters
    **module_name** (`str`)
    : The name of the module to import.

    **ignore_import_errors** (`bool`, optional)
    : Whether to ignore import errors. Default is False.

    #### Yields
    **module** (module)
    : The imported module in the package tree.

    #### Notes
    This function imports the given module and yields it. If `ignore_import_errors`
    is set to True, the function will continue executing even if the import fails
    and will print the traceback. If `ignore_import_errors` is set to False and
    the import fails, the function will raise the error. After yielding the
    imported module, the function looks for sub-modules within the package of
    the imported module and recursively imports and yields them.
    """
    if not ignore_import_errors:
        module = importlib.import_module(module_name)
    else:
        try:
            module = importlib.import_module(module_name)
        except BaseException:
            traceback.print_exc()
            return
    yield module

    if getattr(module, "__path__", None):
        for _, name, _ in pkgutil.iter_modules(module.__path__, f"{module_name}."):
            yield from disc_modules(name, ignore_import_errors)


def disc_benchmarks(root, ignore_import_errors=False):
    """
    Discovers all benchmarks in a given directory tree, yielding Benchmark
    objects.

    #### Parameters
    **root** (`str`)
    : The root of the directory tree where the function begins to search for
      benchmarks.

    **ignore_import_errors** (`bool`, optional)
    : Specifies if import errors should be ignored. Default is False.

    #### Yields
    **benchmark** (Benchmark instance or None)
    : A benchmark instance containing the benchmark's name, the function to
      be benchmarked, and its sources if a matching benchmark is found.

    #### Notes
    For each class definition, the function searches for methods with a
    specific name. For each free function, it yields all functions with a
    specific name. The function initially imports all modules and submodules
    in the directory tree using the `disc_modules` function. Then, for each
    imported module, it searches for classes and functions that might be
    benchmarks. If it finds a class, it looks for methods within that class
    that could be benchmarks. If it finds a free function, it considers it as
    a potential benchmark. A potential benchmark is confirmed by the
    `_get_benchmark` function. If this function returns a benchmark instance,
    the instance is yielded.
    """
    root_name = os.path.basename(root)

    for module in disc_modules(root_name, ignore_import_errors=ignore_import_errors):
        for attr_name, module_attr in (
            (k, v) for k, v in module.__dict__.items() if not k.startswith("_")
        ):
            if inspect.isclass(module_attr) and not inspect.isabstract(module_attr):
                for name, class_attr in inspect.getmembers(module_attr):
                    if inspect.isfunction(class_attr) or inspect.ismethod(class_attr):
                        benchmark = _get_benchmark(
                            name, module, module_attr, class_attr
                        )
                        if benchmark is not None:
                            yield benchmark
            elif inspect.isfunction(module_attr):
                benchmark = _get_benchmark(attr_name, module, None, module_attr)
                if benchmark is not None:
                    yield benchmark


def get_benchmark_from_name(root, name, extra_params=None):
    """
    Creates a benchmark from a fully-qualified benchmark name.

    #### Parameters
    **root** (`str`)
    : Path to the root of a benchmark suite.

    **name** (`str`)
    : Fully-qualified name of a specific benchmark.

    **extra_params** (`dict`, optional)
    : Extra parameters to be added to the benchmark.

    #### Returns
    **benchmark** (Benchmark instance)
    : A benchmark instance created from the given fully-qualified benchmark name.

    #### Raises
    **ValueError**
    : If the provided benchmark ID is invalid or if the benchmark could not be found.

    #### Notes
    This function aims to create a benchmark from the given fully-qualified
    name. It splits the name using the "-" character. If "-" is present in the
    name, the string after the "-" is converted to an integer and is considered as
    the parameter index. If "-" is not present, the parameter index is set to
    None.  The function then tries to directly import the benchmark function by
    guessing its import module name. If the benchmark is not found this way, the
    function searches for the benchmark in the directory tree root using
    `disc_benchmarks`. If the benchmark is still not found, it raises a
    ValueError.  If extra parameters are provided, they are added to the
    benchmark.
    """
    if "-" in name:
        try:
            name, param_idx = name.split("-", 1)
            param_idx = int(param_idx)
        except ValueError:
            raise ValueError(f"Benchmark id {name!r} is invalid")
    else:
        param_idx = None

    update_sys_path(root)
    benchmark = None

    # try to directly import benchmark function by guessing its import module name
    parts = name.split(".")
    for i in [1, 2]:
        path = f"{os.path.join(root, *parts[:-i])}.py"
        if not os.path.isfile(path):
            continue
        modname = ".".join([os.path.basename(root)] + parts[:-i])
        module = importlib.import_module(modname)
        try:
            module_attr = getattr(module, parts[-i])
        except AttributeError:
            break
        if i == 1 and inspect.isfunction(module_attr):
            benchmark = _get_benchmark(parts[-i], module, None, module_attr)
            break
        elif i == 2 and inspect.isclass(module_attr):
            try:
                class_attr = getattr(module_attr, parts[-1])
            except AttributeError:
                break
            if inspect.isfunction(class_attr) or inspect.ismethod(class_attr):
                benchmark = _get_benchmark(parts[-1], module, module_attr, class_attr)
                break

    if benchmark is None:
        for benchmark in disc_benchmarks(root):
            if benchmark.name == name:
                break
        else:
            raise ValueError(f"Could not find benchmark '{name}'")

    if param_idx is not None:
        benchmark.set_param_idx(param_idx)

    if extra_params:

        class ExtraBenchmarkAttrs:
            pass

        for key, value in extra_params.items():
            setattr(ExtraBenchmarkAttrs, key, value)
        benchmark._attr_sources.insert(0, ExtraBenchmarkAttrs)

    return benchmark


def list_benchmarks(root, fp):
    """
    Lists all discovered benchmarks to a file pointer as JSON.

    #### Parameters
    **root** (`str`)
    : Path to the root of a benchmark suite.

    **fp** (file object)
    : File pointer where the JSON list of benchmarks should be written.

    #### Notes
    The function updates the system path with the root directory of the
    benchmark suite. Then, it iterates over all benchmarks discovered in the
    root directory. For each benchmark, it creates a dictionary containing all
    attributes of the benchmark that are of types `str`, `int`, `float`, `list`,
    `dict`, `bool` and don't start with an underscore `_`.  These attribute
    dictionaries are then dumped as JSON into the file pointed by `fp`.
    """
    update_sys_path(root)

    # Streaming of JSON back out to the master process
    fp.write("[")
    first = True
    for benchmark in disc_benchmarks(root):
        if not first:
            fp.write(", ")
        clean = {
            k: v
            for (k, v) in benchmark.__dict__.items()
            if isinstance(v, (str, int, float, list, dict, bool))
            and not k.startswith("_")
        }
        json.dump(clean, fp, skipkeys=True)
        first = False
    fp.write("]")


def _discover(args):
    """
    Discovers all benchmarks in the provided benchmark directory and lists them
    to a file.

    #### Parameters
    **args** (`tuple`)
    : A tuple containing benchmark directory and result file path.

    #### Notes
    The function takes a tuple as an argument. The first element of the tuple
    should be the path to the benchmark directory, and the second element should
    be the path to the result file. It opens the result file for writing and
    calls the `list_benchmarks` function with the benchmark directory and the
    file pointer of the result file.
    """
    benchmark_dir, result_file = args
    with open(result_file, "w") as fp:
        list_benchmarks(benchmark_dir, fp)
