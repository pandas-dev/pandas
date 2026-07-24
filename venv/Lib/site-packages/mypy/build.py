"""Facilities to analyze entire programs, including imported modules.

Parse and analyze the source files of a program in the correct order
(based on file dependencies), and collect the results.

This module only directs a build, which is performed in multiple passes per
file.  The individual passes are implemented in separate modules.

The function build() is the main interface to this module.
"""

# TODO: More consistent terminology, e.g. path/fnam, module/id, state/file

from __future__ import annotations

import collections
import contextlib
import gc
import json
import os
import platform
import re
import stat
import subprocess
import sys
import time
import types
from collections.abc import Callable, Iterator, Mapping, Sequence, Set as AbstractSet
from concurrent.futures import ThreadPoolExecutor, wait
from heapq import heappop, heappush
from textwrap import dedent
from threading import Lock, Thread
from typing import (
    TYPE_CHECKING,
    Any,
    ClassVar,
    Final,
    NoReturn,
    TextIO,
    TypeAlias as _TypeAlias,
    TypedDict,
    cast,
    final,
)

from librt.internal import (
    cache_version,
    read_bool,
    read_int as read_int_bare,
    read_str as read_str_bare,
    read_tag,
    write_bool,
    write_bytes as write_bytes_bare,
    write_int as write_int_bare,
    write_str as write_str_bare,
    write_tag,
)

import mypy.semanal_main
from mypy.cache import (
    CACHE_VERSION,
    DICT_STR_GEN,
    LIST_GEN,
    LITERAL_NONE,
    CacheMeta,
    CacheMetaEx,
    ErrorTuple,
    JsonValue,
    ReadBuffer,
    Tag,
    WriteBuffer,
    read_bytes,
    read_int,
    read_int_list,
    read_str,
    read_str_list,
    read_str_opt,
    write_bytes,
    write_int,
    write_int_list,
    write_json_value,
    write_str,
    write_str_list,
    write_str_opt,
)
from mypy.checker import DeferredNode, TypeChecker
from mypy.defaults import (
    WORKER_CONNECTION_TIMEOUT,
    WORKER_DONE_TIMEOUT,
    WORKER_SHUTDOWN_TIMEOUT,
    WORKER_START_INTERVAL,
    WORKER_START_TIMEOUT,
)
from mypy.error_formatter import OUTPUT_CHOICES, ErrorFormatter
from mypy.errorcodes import ErrorCode
from mypy.errors import CompileError, ErrorInfo, Errors, report_internal_error
from mypy.graph_utils import prepare_sccs, strongly_connected_components, topsort
from mypy.indirection import TypeIndirectionVisitor
from mypy.ipc import (
    BadStatus,
    IPCClient,
    IPCException,
    IPCMessage,
    read_status,
    ready_to_read,
    receive,
    send,
)
from mypy.messages import MessageBuilder
from mypy.nodes import (
    Decorator,
    FileRawData,
    FuncDef,
    Import,
    ImportAll,
    ImportBase,
    ImportFrom,
    MypyFile,
    OverloadedFuncDef,
    SymbolTable,
)
from mypy.options import OPTIONS_AFFECTING_CACHE_NO_PLATFORM
from mypy.partially_defined import PossiblyUndefinedVariableVisitor
from mypy.semanal import SemanticAnalyzer
from mypy.semanal_pass1 import SemanticAnalyzerPreAnalysis
from mypy.util import (
    DecodeError,
    decode_python_encoding,
    get_available_threads,
    get_mypy_comments,
    hash_digest,
    hash_digest_bytes,
    is_stub_package_file,
    is_sub_path_normabs,
    is_typeshed_file,
    module_prefix,
    os_path_join,
    read_py_file,
    time_ref,
    time_spent_us,
)

if TYPE_CHECKING:
    from mypy.report import Reports  # Avoid unconditional slow import

from mypy import errorcodes as codes
from mypy.config_parser import get_config_module_names, parse_mypy_comments
from mypy.fixup import NodeFixer
from mypy.freetree import free_tree
from mypy.fscache import FileSystemCache
from mypy.known_modules import get_known_modules, reset_known_modules_cache
from mypy.messages import best_matches, pretty_seq
from mypy.metastore import FilesystemMetadataStore, MetadataStore, SqliteMetadataStore
from mypy.modulefinder import (
    BuildSource as BuildSource,
    BuildSourceSet as BuildSourceSet,
    FindModuleCache,
    ModuleNotFoundReason,
    ModuleSearchResult,
    SearchPaths,
    compute_search_paths,
)
from mypy.modules_state import modules_state
from mypy.nodes import Expression
from mypy.options import Options
from mypy.parse import load_from_raw, parse
from mypy.plugin import ChainedPlugin, Plugin, ReportConfigContext
from mypy.plugins.default import DefaultPlugin
from mypy.renaming import LimitedVariableRenameVisitor, VariableRenameVisitor
from mypy.stats import dump_type_stats
from mypy.stubinfo import stub_distribution_name
from mypy.types import Type, instance_cache
from mypy.typestate import reset_global_state, type_state
from mypy.util import json_dumps, json_loads
from mypy.version import __version__

# Switch to True to produce debug output related to fine-grained incremental
# mode only that is useful during development. This produces only a subset of
# output compared to --verbose output. We use a global flag to enable this so
# that it's easy to enable this when running tests.
DEBUG_FINE_GRAINED: Final = False

# These modules are special and should always come from typeshed.
CORE_BUILTIN_MODULES: Final = {
    "builtins",
    "typing",
    "types",
    "typing_extensions",
    "mypy_extensions",
    "_typeshed",
    "_collections_abc",
    "collections",
    "collections.abc",
    "sys",
    "abc",
}

# We are careful now, we can increase this in future if safe/useful.
MAX_GC_FREEZE_CYCLES: Final = 1

# We store status of initial GC freeze as a global variable to avoid memory
# leaks in tests, where we keep creating new BuildManagers in the same process.
initial_gc_freeze_done = False

Graph: _TypeAlias = dict[str, "State"]

MODULE_RESOLUTION_URL: Final = (
    "https://mypy.readthedocs.io/en/stable/running_mypy.html#mapping-file-paths-to-modules"
)

# Padding when estimating how much time it will take to process a file. This is to avoid
# situations where 100 empty __init__.py files cost less than 1 trivial module.
MIN_SIZE_HINT: Final = 256


class SCC:
    """A simple class that represents a strongly connected component (import cycle)."""

    id_counter: ClassVar[int] = 0

    def __init__(
        self, ids: set[str], scc_id: int | None = None, deps: list[int] | None = None
    ) -> None:
        if scc_id is None:
            self.id = SCC.id_counter
            SCC.id_counter += 1
        else:
            self.id = scc_id
        # Ids of modules in this cycle.
        self.mod_ids = ids
        # Direct dependencies, should be populated by the caller.
        self.deps: set[int] = set(deps) if deps is not None else set()
        # Direct dependencies that have not been processed yet.
        # Should be populated by the caller. This set may change during graph
        # processing, while the above stays constant.
        self.not_ready_deps: set[int] = set()
        # SCCs that (directly) depend on this SCC. Note this is a list to
        # make processing order more predictable. Dependents will be notified
        # that they may be ready in the order in this list.
        self.direct_dependents: list[int] = []
        # Rough estimate of how much time processing this SCC will take, this
        # is used for more efficient scheduling across multiple build workers.
        self.size_hint: int = MIN_SIZE_HINT


# TODO: Get rid of BuildResult.  We might as well return a BuildManager.
class BuildResult:
    """The result of a successful build.

    Attributes:
      manager: The build manager.
      files:   Dictionary from module name to related AST node.
      types:   Dictionary from parse tree node to its inferred type.
      used_cache: Whether the build took advantage of a pre-existing cache
      errors:  List of error messages.
    """

    def __init__(self, manager: BuildManager, graph: Graph) -> None:
        self.manager = manager
        self.graph = graph
        self.files = manager.modules
        self.types = manager.all_types  # Non-empty if export_types True in options
        self.used_cache = manager.cache_enabled
        self.errors: list[str] = []  # Filled in by build if desired


class WorkerClient:
    """A simple class that represents a mypy build worker."""

    conn: IPCClient

    def __init__(self, status_file: str, options_data: str, env: Mapping[str, str]) -> None:
        self.status_file = status_file
        if os.path.isfile(status_file):
            os.unlink(status_file)

        command = [
            sys.executable,
            "-m",
            "mypy.build_worker",
            f"--status-file={status_file}",
            f"--options-data={options_data}",
        ]
        # Return early without waiting, caller must call connect() before using the client.
        self.proc = subprocess.Popen(command, env=env)
        self.connected = False

    def connect(self) -> None:
        end_time = time.time() + WORKER_START_TIMEOUT
        last_exception: Exception | None = None
        while time.time() < end_time:
            try:
                data = read_status(self.status_file)
            except BadStatus as exc:
                last_exception = exc
                time.sleep(WORKER_START_INTERVAL)
                continue
            try:
                pid, connection_name = data["pid"], data["connection_name"]
                assert isinstance(pid, int), f"Bad PID: {pid}"
                assert isinstance(connection_name, str), f"Bad connection name: {connection_name}"
                if sys.platform != "win32":
                    # Windows uses "wrapper processes" to run Python, so we cannot
                    # verify PIDs reliably.
                    assert pid == self.proc.pid, f"PID mismatch: {pid} vs {self.proc.pid}"
                self.conn = IPCClient(connection_name, WORKER_CONNECTION_TIMEOUT)
                self.connected = True
                return
            except Exception as exc:
                last_exception = exc
                break
        print(f"Failed to establish connection with worker: {last_exception}")

    def close(self) -> None:
        if self.connected:
            self.conn.close()
        # Technically we don't need to wait, but otherwise we will get ResourceWarnings.
        try:
            self.proc.wait(timeout=WORKER_SHUTDOWN_TIMEOUT)
        except subprocess.TimeoutExpired:
            pass
        if os.path.isfile(self.status_file):
            os.unlink(self.status_file)


def build_error(msg: str) -> NoReturn:
    raise CompileError([f"mypy: error: {msg}"])


def build(
    sources: list[BuildSource],
    options: Options,
    alt_lib_path: str | None = None,
    flush_errors: Callable[[str | None, list[str], bool], None] | None = None,
    fscache: FileSystemCache | None = None,
    stdout: TextIO | None = None,
    stderr: TextIO | None = None,
    extra_plugins: Sequence[Plugin] | None = None,
    worker_env: Mapping[str, str] | None = None,
) -> BuildResult:
    """Analyze a program.

    A single call to build performs parsing, semantic analysis and optionally
    type checking for the program *and* all imported modules, recursively.

    Return BuildResult if successful or only non-blocking errors were found;
    otherwise raise CompileError.

    If a flush_errors callback is provided, all error messages will be
    passed to it and the errors and messages fields of BuildResult and
    CompileError (respectively) will be empty. Otherwise, those fields will
    report any error messages.

    Args:
      sources: list of sources to build
      options: build options
      alt_lib_path: an additional directory for looking up library modules
        (takes precedence over other directories)
      flush_errors: optional function to flush errors after a file is processed
      fscache: optionally a file-system cacher
      stdout: Output stream to use instead of `sys.stdout`
      stderr: Error stream to use instead of `sys.stderr`
      extra_plugins: Plugins to use in addition to those loaded from config
      worker_env: An environment to start parallel build workers (used for tests)
    """
    # If we were not given a flush_errors, we use one that will populate those
    # fields for callers that want the traditional API.
    messages = []

    # This is mostly for the benefit of tests that use builtins fixtures.
    instance_cache.reset()
    reset_known_modules_cache()

    def default_flush_errors(
        filename: str | None, new_messages: list[str], is_serious: bool
    ) -> None:
        messages.extend(new_messages)

    flush_errors = flush_errors or default_flush_errors
    stdout = stdout or sys.stdout
    stderr = stderr or sys.stderr
    extra_plugins = extra_plugins or []

    # Create metastore before workers to avoid race conditions.
    metastore = create_metastore(options, parallel_worker=False)
    workers = []
    connect_threads = []
    # A quasi-unique ID for this specific mypy invocation.
    build_id = os.urandom(4).hex()
    options_data = None
    if options.num_workers > 0:
        os.makedirs(options.cache_dir, exist_ok=True)
        options_data = os_path_join(options.cache_dir, f".worker_options.{build_id}.data")
        with open(options_data, "wb") as f:
            f.write(options.to_bytes())
        workers = [
            WorkerClient(
                f".mypy_worker.{build_id}.{idx}.json", options_data, worker_env or os.environ
            )
            for idx in range(options.num_workers)
        ]
        sources_message = SourcesDataMessage(sources=sources)
        buf = WriteBuffer()
        sources_message.write(buf)
        sources_data = buf.getvalue()

        def connect(wc: WorkerClient, data: bytes) -> None:
            # Start loading sources in each worker as soon as it is up.
            wc.connect()
            if not wc.connected:
                # Caller should detect this and fail gracefully.
                return
            wc.conn.write_bytes(data)

        # We don't wait for workers to be ready until they are actually needed.
        for worker in workers:
            thread = Thread(target=connect, args=(worker, sources_data))
            thread.start()
            connect_threads.append(thread)

    try:
        result = build_inner(
            sources,
            options,
            alt_lib_path,
            flush_errors,
            fscache,
            stdout,
            stderr,
            extra_plugins,
            workers,
            connect_threads,
            metastore,
        )
        result.errors = messages
        return result
    except CompileError as e:
        # CompileErrors raised from an errors object carry all the
        # messages that have not been reported out by error streaming.
        # Patch it up to contain either none or all none of the messages,
        # depending on whether we are flushing errors.
        serious = not e.use_stdout
        flush_errors(None, e.messages, serious)
        e.messages = messages
        raise
    finally:
        # In case of an early crash it is better to wait for workers to become ready, and
        # shut them down cleanly. Otherwise, they will linger until connection timeout.
        for thread in connect_threads:
            thread.join()
        if options_data is not None:
            os.unlink(options_data)
        for worker in workers:
            if not worker.connected:
                continue
            try:
                send(worker.conn, SccRequestMessage(scc_ids=[], import_errors={}, mod_data={}))
            except (OSError, IPCException):
                pass
        for worker in workers:
            worker.close()


def build_inner(
    sources: list[BuildSource],
    options: Options,
    alt_lib_path: str | None,
    flush_errors: Callable[[str | None, list[str], bool], None],
    fscache: FileSystemCache | None,
    stdout: TextIO,
    stderr: TextIO,
    extra_plugins: Sequence[Plugin],
    workers: list[WorkerClient],
    connect_threads: list[Thread],
    metastore: MetadataStore,
) -> BuildResult:
    if platform.python_implementation() == "CPython":
        # Run gc less frequently, as otherwise we can spend a large fraction of
        # cpu in gc. This seems the most reasonable place to tune garbage collection.
        gc.set_threshold(200 * 1000, 30, 30)

    data_dir = default_data_dir()
    fscache = fscache or FileSystemCache()

    search_paths = compute_search_paths(sources, options, data_dir, alt_lib_path)

    reports = None
    if options.report_dirs:
        # Import lazily to avoid slowing down startup.
        from mypy.report import Reports

        reports = Reports(data_dir, options.report_dirs)

    source_set = BuildSourceSet(sources)
    cached_read = fscache.read
    error_formatter = None if options.output is None else OUTPUT_CHOICES.get(options.output)
    errors = Errors(
        options,
        read_source=lambda path: read_py_file(path, cached_read),
        error_formatter=error_formatter,
    )
    # Record import errors so that they can be replayed by the workers.
    if workers:
        errors.global_watcher = True
    plugin, snapshot = load_plugins(options, errors, stdout, extra_plugins)

    # Validate error codes after plugins are loaded.
    options.process_error_codes(error_callback=build_error)

    # Construct a build manager object to hold state during the build.
    #
    # Ignore current directory prefix in error messages.
    manager = BuildManager(
        data_dir,
        search_paths,
        ignore_prefix=os.getcwd(),
        source_set=source_set,
        reports=reports,
        options=options,
        version_id=__version__,
        plugin=plugin,
        plugins_snapshot=snapshot,
        errors=errors,
        error_formatter=error_formatter,
        flush_errors=flush_errors,
        fscache=fscache,
        stdout=stdout,
        stderr=stderr,
        metastore=metastore,
    )
    manager.workers = workers
    if manager.verbosity() >= 2:
        manager.trace(repr(options))

    reset_global_state()
    try:
        graph = dispatch(sources, manager, stdout, connect_threads)
        if not options.fine_grained_incremental:
            type_state.reset_all_subtype_caches()
        if options.timing_stats is not None:
            dump_timing_stats(options.timing_stats, graph)
        if options.line_checking_stats is not None:
            dump_line_checking_stats(options.line_checking_stats, graph)
        warn_unused_configs(options, flush_errors)
        return BuildResult(manager, graph)
    finally:
        manager.commit()
        manager.log(
            "Build finished in %.3f seconds with %d modules, and %d errors"
            % (
                time.time() - manager.start_time,
                len(manager.modules),
                manager.errors.num_messages(),
            )
        )
        manager.dump_stats()
        if reports is not None:
            # Finish the HTML or XML reports even if CompileError was raised.
            reports.finish()
        if os.path.isdir(options.cache_dir):
            add_catch_all_gitignore(options.cache_dir)
            exclude_from_backups(options.cache_dir)
        if os.path.isdir(options.cache_dir):
            record_missing_stub_packages(options.cache_dir, manager.missing_stub_packages)


def warn_unused_configs(
    options: Options, flush_errors: Callable[[str | None, list[str], bool], None]
) -> None:
    unused_configs = options.get_unused_configs()
    if options.warn_unused_configs and unused_configs and not options.non_interactive:
        unused = get_config_module_names(
            options.config_file,
            [glob for glob in options.per_module_options.keys() if glob in unused_configs],
        )
        flush_errors(
            None, ["{}: note: unused section(s): {}".format(options.config_file, unused)], False
        )


def default_data_dir() -> str:
    """Returns directory containing typeshed directory."""
    return os.path.dirname(__file__)


def normpath(path: str, options: Options) -> str:
    """Convert path to absolute; but to relative in bazel mode.

    (Bazel's distributed cache doesn't like filesystem metadata to
    end up in output files.)
    """
    # TODO: Could we always use relpath?  (A worry in non-bazel
    # mode would be that a moved file may change its full module
    # name without changing its size, mtime or hash.)
    if options.bazel:
        return os.path.relpath(path)
    else:
        return os.path.abspath(path)


# NOTE: dependencies + suppressed == all reachable imports;
# suppressed contains those reachable imports that were prevented by
# silent mode or simply not found.


# Metadata for the fine-grained dependencies file associated with a module.
class FgDepMeta(TypedDict):
    path: str
    mtime: int


# Priorities used for imports.  (Here, top-level includes inside a class.)
# These are used to determine a more predictable order in which the
# nodes in an import cycle are processed.
PRI_HIGH: Final = 5  # top-level "from X import blah"
PRI_MED: Final = 10  # top-level "import X"
PRI_LOW: Final = 20  # either form inside a function
PRI_MYPY: Final = 25  # inside "if MYPY" or "if TYPE_CHECKING"
PRI_INDIRECT: Final = 30  # an indirect dependency
PRI_ALL: Final = 99  # include all priorities


def import_priority(imp: ImportBase, toplevel_priority: int) -> int:
    """Compute import priority from an import node."""
    if not imp.is_top_level:
        # Inside a function
        return PRI_LOW
    if imp.is_mypy_only:
        # Inside "if MYPY" or "if typing.TYPE_CHECKING"
        return max(PRI_MYPY, toplevel_priority)
    # A regular import; priority determined by argument.
    return toplevel_priority


def load_plugins_from_config(
    options: Options, errors: Errors, stdout: TextIO
) -> tuple[list[Plugin], dict[str, str]]:
    """Load all configured plugins.

    Return a list of all the loaded plugins from the config file.
    The second return value is a snapshot of versions/hashes of loaded user
    plugins (for cache validation).
    """
    import importlib

    snapshot: dict[str, str] = {}

    if not options.config_file:
        return [], snapshot

    line = find_config_file_line_number(options.config_file, "mypy", "plugins")
    if line == -1:
        line = 1  # We need to pick some line number that doesn't look too confusing

    def plugin_error(message: str) -> NoReturn:
        errors.report(line, 0, message)
        errors.raise_error(use_stdout=False)

    custom_plugins: list[Plugin] = []
    errors.set_file(options.config_file, None, options)
    for plugin_path in options.plugins:
        func_name = "plugin"
        plugin_dir: str | None = None
        if ":" in os.path.basename(plugin_path):
            plugin_path, func_name = plugin_path.rsplit(":", 1)
        if plugin_path.endswith(".py"):
            # Plugin paths can be relative to the config file location.
            plugin_path = os_path_join(os.path.dirname(options.config_file), plugin_path)
            if not os.path.isfile(plugin_path):
                plugin_error(f'Can\'t find plugin "{plugin_path}"')
            # Use an absolute path to avoid populating the cache entry
            # for 'tmp' during tests, since it will be different in
            # different tests.
            plugin_dir = os.path.abspath(os.path.dirname(plugin_path))
            fnam = os.path.basename(plugin_path)
            module_name = fnam[:-3]
            sys.path.insert(0, plugin_dir)
        elif re.search(r"[\\/]", plugin_path):
            fnam = os.path.basename(plugin_path)
            plugin_error(f'Plugin "{fnam}" does not have a .py extension')
        else:
            module_name = plugin_path

        try:
            module = importlib.import_module(module_name)
        except Exception as exc:
            plugin_error(f'Error importing plugin "{plugin_path}": {exc}')
        finally:
            if plugin_dir is not None:
                assert sys.path[0] == plugin_dir
                del sys.path[0]

        if not hasattr(module, func_name):
            plugin_error(
                'Plugin "{}" does not define entry point function "{}"'.format(
                    plugin_path, func_name
                )
            )

        try:
            plugin_type = getattr(module, func_name)(__version__)
        except Exception:
            print(f"Error calling the plugin(version) entry point of {plugin_path}\n", file=stdout)
            raise  # Propagate to display traceback

        if not isinstance(plugin_type, type):
            plugin_error(
                'Type object expected as the return value of "plugin"; got {!r} (in {})'.format(
                    plugin_type, plugin_path
                )
            )
        if not issubclass(plugin_type, Plugin):
            plugin_error(
                'Return value of "plugin" must be a subclass of "mypy.plugin.Plugin" '
                "(in {})".format(plugin_path)
            )
        try:
            custom_plugins.append(plugin_type(options))
            snapshot[module_name] = take_module_snapshot(module)
        except Exception:
            print(f"Error constructing plugin instance of {plugin_type.__name__}\n", file=stdout)
            raise  # Propagate to display traceback

    return custom_plugins, snapshot


def load_plugins(
    options: Options, errors: Errors, stdout: TextIO, extra_plugins: Sequence[Plugin]
) -> tuple[Plugin, dict[str, str]]:
    """Load all configured plugins.

    Return a plugin that encapsulates all plugins chained together. Always
    at least include the default plugin (it's last in the chain).
    The second return value is a snapshot of versions/hashes of loaded user
    plugins (for cache validation).
    """
    custom_plugins, snapshot = load_plugins_from_config(options, errors, stdout)

    custom_plugins += extra_plugins

    default_plugin: Plugin = DefaultPlugin(options)
    if not custom_plugins:
        return default_plugin, snapshot

    # Custom plugins take precedence over the default plugin.
    return ChainedPlugin(options, custom_plugins + [default_plugin]), snapshot


def take_module_snapshot(module: types.ModuleType) -> str:
    """Take plugin module snapshot by recording its version and hash.

    We record _both_ hash and the version to detect more possible changes
    (e.g. if there is a change in modules imported by a plugin).
    """
    if hasattr(module, "__file__"):
        assert module.__file__ is not None
        with open(module.__file__, "rb") as f:
            digest = hash_digest(f.read())
    else:
        digest = "unknown"
    ver = getattr(module, "__version__", "none")
    return f"{ver}:{digest}"


def find_config_file_line_number(path: str, section: str, setting_name: str) -> int:
    """Return the approximate location of setting_name within mypy config file.

    Return -1 if can't determine the line unambiguously.
    """
    in_desired_section = False
    try:
        results = []
        with open(path, encoding="UTF-8") as f:
            for i, line in enumerate(f):
                line = line.strip()
                if line.startswith("[") and line.endswith("]"):
                    current_section = line[1:-1].strip()
                    in_desired_section = current_section == section
                elif in_desired_section and re.match(rf"{setting_name}\s*=", line):
                    results.append(i + 1)
        if len(results) == 1:
            return results[0]
    except OSError:
        pass
    return -1


class BuildManager:
    """This class holds shared state for building a mypy program.

    It is used to coordinate parsing, import processing, semantic
    analysis and type checking.  The actual build steps are carried
    out by dispatch().

    Attributes:
      data_dir:        Mypy data directory (contains stubs)
      search_paths:    SearchPaths instance indicating where to look for modules
      modules:         Mapping of module ID to MypyFile (shared by the passes)
      semantic_analyzer:
                       Semantic analyzer, pass 2
      all_types:       Map {Expression: Type} from all modules (enabled by export_types)
      options:         Build options
      missing_modules: Modules that could not be imported (or intentionally skipped)
      stale_modules:   Set of modules that needed to be rechecked (only used by tests)
      fg_deps_meta:    Metadata for fine-grained dependencies caches associated with modules
      fg_deps:         A fine-grained dependency map
      version_id:      The current mypy version (based on commit id when possible)
      plugin:          Active mypy plugin(s)
      plugins_snapshot:
                       Snapshot of currently active user plugins (versions and hashes)
      old_plugins_snapshot:
                       Plugins snapshot from previous incremental run (or None in
                       non-incremental mode and if cache was not found)
      errors:          Used for reporting all errors
      flush_errors:    A function for processing errors after each SCC
      cache_enabled:   Whether cache is being read. This is set based on options,
                       but is disabled if fine-grained cache loading fails
                       and after an initial fine-grained load. This doesn't
                       determine whether we write cache files or not.
      quickstart_state:
                       A cache of filename -> mtime/size/hash info used to avoid
                       needing to hash source files when using a cache with mismatching mtimes
      stats:           Dict with various instrumentation numbers, it is used
                       not only for debugging, but also required for correctness,
                       in particular to check consistency of the fine-grained dependency cache.
      fscache:         A file system cacher
      ast_cache:       AST cache to speed up mypy daemon
    """

    def __init__(
        self,
        data_dir: str,
        search_paths: SearchPaths,
        ignore_prefix: str,
        source_set: BuildSourceSet,
        reports: Reports | None,
        options: Options,
        version_id: str,
        plugin: Plugin,
        plugins_snapshot: dict[str, str],
        errors: Errors,
        flush_errors: Callable[[str | None, list[str], bool], None],
        fscache: FileSystemCache,
        stdout: TextIO,
        stderr: TextIO,
        error_formatter: ErrorFormatter | None = None,
        parallel_worker: bool = False,
        metastore: MetadataStore | None = None,
    ) -> None:
        self.stats: dict[str, Any] = {}  # Values are ints or floats
        # Use in cases where we need to prevent race conditions in stats reporting.
        self.stats_lock = Lock()
        self.stdout = stdout
        self.stderr = stderr
        self.start_time = time.time()
        self.data_dir = data_dir
        self.errors = errors
        self.errors.set_ignore_prefix(ignore_prefix)
        self.error_formatter = error_formatter
        self.search_paths = search_paths
        self.source_set = source_set
        self.reports = reports
        self.options = options
        self.version_id = version_id
        self.modules: dict[str, MypyFile] = {}
        # Share same modules dictionary with the global fixer state.
        # We need to set allow_missing when doing a fine-grained cache
        # load because we need to gracefully handle missing modules.
        modules_state.modules = self.modules
        modules_state.node_fixer = NodeFixer(self.modules, self.options.use_fine_grained_cache)
        self.import_map: dict[str, set[str]] = {}
        self.missing_modules: dict[str, int] = {}
        self.fg_deps_meta: dict[str, FgDepMeta] = {}
        # fg_deps holds the dependencies of every module that has been
        # processed. We store this in BuildManager so that we can compute
        # dependencies as we go, which allows us to free ASTs and type information,
        # saving a ton of memory on net.
        self.fg_deps: dict[str, set[str]] = {}
        # Always convert the plugin to a ChainedPlugin so that it can be manipulated if needed
        if not isinstance(plugin, ChainedPlugin):
            plugin = ChainedPlugin(options, [plugin])
        self.plugin = plugin
        # These allow quickly skipping logging and stats collection calls. Note
        # that some stats impact mypy behavior, so be careful when skipping stats
        # collection calls.
        self.stats_enabled = self.options.dump_build_stats
        self.logging_enabled = self.options.verbosity >= 1
        self.tracing_enabled = self.options.verbosity >= 2
        # Set of namespaces (module or class) that are being populated during semantic
        # analysis and may have missing definitions.
        self.incomplete_namespaces: set[str] = set()
        self.semantic_analyzer = SemanticAnalyzer(
            self.modules,
            self.missing_modules,
            self.incomplete_namespaces,
            self.errors,
            self.plugin,
            self.import_map,
        )
        self.all_types: dict[Expression, Type] = {}  # Enabled by export_types
        self.indirection_detector = TypeIndirectionVisitor()
        self.stale_modules: set[str] = set()
        self.rechecked_modules: set[str] = set()
        self.flush_errors = flush_errors
        has_reporters = reports is not None and reports.reporters
        self.cache_enabled = (
            options.incremental
            and (not options.fine_grained_incremental or options.use_fine_grained_cache)
            and not has_reporters
        )
        self.fscache = fscache
        self.cwd = os.getcwd()
        self.find_module_cache = FindModuleCache(
            self.search_paths, self.fscache, self.options, source_set=self.source_set
        )
        for module in CORE_BUILTIN_MODULES:
            if options.use_builtins_fixtures:
                continue
            path = self.find_module_cache.find_module(module, fast_path=True)
            if not isinstance(path, str):
                build_error(
                    f'Failed to find builtin module "{module}", perhaps typeshed is broken?'
                )
            if is_typeshed_file(options.abs_custom_typeshed_dir, path) or is_stub_package_file(
                path
            ):
                continue

            self.errors.set_file(path, module, options)
            self.error(None, f'This file shadows library module "{module}"', blocker=True)
            self.note(
                None, f'A user-defined top-level module with name "{module}" is not supported'
            )
            self.errors.raise_error()

        if metastore is None:
            metastore = create_metastore(options, parallel_worker=parallel_worker)
        self.metastore = metastore

        # a mapping from source files to their corresponding shadow files
        # for efficient lookup
        self.shadow_map: dict[str, str] = {}
        if self.options.shadow_file is not None:
            self.shadow_map = dict(self.options.shadow_file)
        # a mapping from each file being typechecked to its possible shadow file
        self.shadow_equivalence_map: dict[str, str | None] = {}
        self.plugin = plugin
        self.plugins_snapshot = plugins_snapshot
        self.old_plugins_snapshot = read_plugins_snapshot(self)
        if self.verbosity() >= 2:
            self.trace(f"Plugins snapshot (fresh) {json.dumps(self.plugins_snapshot)}")
        self.quickstart_state = read_quickstart_file(options, self.stdout)
        # Fine grained targets (module top levels and top level functions) processed by
        # the semantic analyzer, used only for testing. Currently used only by the new
        # semantic analyzer. Tuple of module and target name.
        self.processed_targets: list[tuple[str, str]] = []
        # Missing stub packages encountered.
        self.missing_stub_packages: set[str] = set()
        # Cache for mypy ASTs that have completed semantic analysis
        # pass 1. When multiple files are added to the build in a
        # single daemon increment, only one of the files gets added
        # per step and the others are discarded. This gets repeated
        # until all the files have been added. This means that a
        # new file can be processed O(n**2) times. This cache
        # avoids most of this redundant work.
        self.ast_cache: dict[str, tuple[MypyFile, list[ErrorInfo], str | None]] = {}
        # Number of times we used GC optimization hack for fresh SCCs.
        self.gc_freeze_cycles = 0
        # Mapping from SCC id to corresponding SCC instance. This is populated
        # in process_graph().
        self.scc_by_id: dict[int, SCC] = {}
        # Mapping from module id to the SCC it belongs to. This is populated
        # in process_graph().
        self.scc_by_mod_id: dict[str, SCC] = {}
        # Global topological order for SCCs. This exists to make order of processing
        # SCCs more predictable.
        self.top_order: list[int] = []
        # Stale SCCs that are queued for processing. Each tuple contains SCC size hint,
        # SCC adding order (tie-breaker), and the SCC itself.
        self.scc_queue: list[tuple[int, int, SCC]] = []
        # Total size hint for SCCs currently in queue.
        self.size_in_queue: int = 0
        # SCCs that have been fully processed.
        self.done_sccs: set[int] = set()
        # Parallel build workers, list is empty for in-process type-checking.
        self.workers: list[WorkerClient] = []
        # We track which workers are currently free in the coordinator process.
        # This is a tiny bit faster and conceptually simpler than check which ones
        # are writeable each time we want to submit an SCC for processing.
        self.free_workers: set[int] = set()
        # A global adding order for SCC queue, see comment above.
        self.queue_order: int = 0
        # Is this an instance used by a parallel worker?
        self.parallel_worker = parallel_worker
        # Snapshot of import-related options per module. We record these even for
        # suppressed imports, since they can affect errors in the callers. Bytes
        # value is opaque but can be compared to detect changes in options.
        self.import_options: dict[str, bytes] = {}
        # Cache for transitive dependency check (expensive).
        self.transitive_deps_cache: dict[tuple[int, int], bool] = {}
        # Packages for which we know presence or absence of __getattr__().
        self.known_partial_packages: dict[str, bool] = {}

    def dump_stats(self) -> None:
        if self.stats_enabled:
            lines = ["Stats:"]
            for key, value in sorted(self.stats_summary().items()):
                fmt = ".3f" if isinstance(value, float) else "d"
                lines.append(f"{key + ':':24}{value:{fmt}}")
            # Call print once so that we don't get a mess in parallel mode.
            print("\n".join(lines) + "\n\n", end="")

    def parse_all(self, states: list[State], post_parse: bool = True) -> None:
        """Parse multiple files in parallel (if possible) and compute dependencies.

        If post_parse is False, skip the last step (used when parsing unchanged files
        that need to be re-checked due to stale dependencies).
        """
        if not self.options.native_parser:
            # Old parser cannot be parallelized.
            for state in states:
                state.parse_file()
            if post_parse:
                self.post_parse_all(states)
            return

        sequential_states = []
        parallel_states = []
        for state in states:
            if state.tree is not None:
                # The file was already parsed.
                continue
            if not self.fscache.exists(state.xpath, real_only=True):
                # New parser only supports parsing on-disk files.
                sequential_states.append(state)
                continue
            parallel_states.append(state)
        if len(parallel_states) > 1:
            self.parse_parallel(sequential_states, parallel_states)
        else:
            # Avoid using executor when there is no parallelism.
            for state in states:
                state.parse_file()
        if post_parse:
            self.post_parse_all(states)

    def parse_parallel(self, sequential_states: list[State], parallel_states: list[State]) -> None:
        """Perform parallel parsing of states.

        Note: this duplicates a bit of logic from State.parse_file(). This is done
        as an optimization to parallelize only those parts of the code that can be
        parallelized efficiently.
        """
        parallel_parsed_states, parallel_parsed_states_set = self.parse_files_threaded_raw(
            sequential_states, parallel_states
        )

        for state in parallel_parsed_states:
            # New parser returns serialized ASTs. Deserialize full trees only if not using
            # parallel workers.
            with state.wrap_context():
                assert state.tree is not None
                raw_data = state.tree.raw_data
                if raw_data is not None:
                    # Apply inline mypy config before deserialization, since
                    # some options (e.g. implicit_optional) affect deserialization
                    state.source_hash = raw_data.source_hash
                    state.apply_inline_configuration(raw_data.mypy_comments)
                    state.tree = load_from_raw(
                        state.xpath,
                        state.id,
                        raw_data,
                        self.errors,
                        state.options,
                        imports_only=bool(self.workers),
                    )
                if self.errors.is_blockers():
                    self.log("Bailing due to parse errors")
                    self.errors.raise_error()

        for state in parallel_states:
            assert state.tree is not None
            if state in parallel_parsed_states_set:
                if state.tree.raw_data is not None:
                    # source_hash was already extracted above, but raw_data
                    # may have been preserved for workers (imports_only=True).
                    pass
                elif state.source_hash is None:
                    # At least namespace packages may not have source.
                    state.get_source()
                state.early_errors = list(self.errors.error_info_map.get(state.xpath, []))
                state.semantic_analysis_pass1()
                self.ast_cache[state.id] = (state.tree, state.early_errors, state.source_hash)
            self.modules[state.id] = state.tree
            if state.tree.raw_data is not None:
                state.size_hint = len(state.tree.raw_data.defs) + MIN_SIZE_HINT
            state.check_blockers()
            state.setup_errors()

    def parse_files_threaded_raw(
        self, sequential_states: list[State], parallel_states: list[State]
    ) -> tuple[list[State], set[State]]:
        """Parse files using a thread pool.

        Also parse sequential states while waiting for the parallel results.
        Trees from the new parser are left in raw (serialized) form.

        Return (list, set) of states that were actually parsed (not cached).
        """
        futures = []
        # Use both list and a set to have more predictable order of errors,
        # while also not sacrificing performance.
        parallel_parsed_states: list[State] = []
        parallel_parsed_states_set: set[State] = set()
        # Use at least --num-workers if specified by user.
        available_threads = max(get_available_threads(), self.options.num_workers)
        # Overhead from trying to parallelize (small) blocking portion of
        # parse_file_inner() results in no visible improvement with more than 8 threads.
        # TODO: reuse thread pool and/or batch small files in single submit() call.
        with ThreadPoolExecutor(max_workers=min(available_threads, 8)) as executor:
            for state in parallel_states:
                state.needs_parse = False
                if state.id not in self.ast_cache:
                    self.log(f"Parsing {state.xpath} ({state.id})")
                    ignore_errors = state.ignore_all or state.options.ignore_errors
                    if ignore_errors:
                        self.errors.ignored_files.add(state.xpath)
                    futures.append(executor.submit(state.parse_file_inner, ""))
                    parallel_parsed_states.append(state)
                    parallel_parsed_states_set.add(state)
                else:
                    self.log(f"Using cached AST for {state.xpath} ({state.id})")
                    state.tree, state.early_errors, source_hash = self.ast_cache[state.id]
                    state.source_hash = source_hash

            # Parse sequential before waiting on parallel.
            for state in sequential_states:
                state.parse_file()

            for fut in wait(futures).done:
                fut.result()

        return parallel_parsed_states, parallel_parsed_states_set

    def post_parse_all(self, states: list[State]) -> None:
        for state in states:
            state.compute_dependencies()
            if self.workers and state.tree:
                # We don't need imports in coordinator process anymore, we parse only to
                # compute dependencies.
                state.tree.imports = []
                del self.ast_cache[state.id]

    def use_fine_grained_cache(self) -> bool:
        return self.cache_enabled and self.options.use_fine_grained_cache

    def maybe_swap_for_shadow_path(self, path: str) -> str:
        if not self.shadow_map:
            return path

        path = normpath(path, self.options)

        previously_checked = path in self.shadow_equivalence_map
        if not previously_checked:
            for source, shadow in self.shadow_map.items():
                if self.fscache.samefile(path, source):
                    self.shadow_equivalence_map[path] = shadow
                    break
                else:
                    self.shadow_equivalence_map[path] = None

        shadow_file = self.shadow_equivalence_map.get(path)
        return shadow_file if shadow_file else path

    def get_stat(self, path: str) -> os.stat_result | None:
        return self.fscache.stat_or_none(self.maybe_swap_for_shadow_path(path))

    def getmtime(self, path: str) -> int:
        """Return a file's mtime; but 0 in bazel mode.

        (Bazel's distributed cache doesn't like filesystem metadata to
        end up in output files.)
        """
        if self.options.bazel:
            return 0
        else:
            return int(self.metastore.getmtime(path))

    def correct_rel_imp(self, file: MypyFile, imp: ImportFrom | ImportAll) -> str:
        """Function to correct for relative imports."""
        file_id = file.fullname
        rel = imp.relative
        if rel == 0:
            return imp.id
        if os.path.basename(file.path).startswith("__init__."):
            rel -= 1
        if rel != 0:
            file_id = ".".join(file_id.split(".")[:-rel])
        new_id = file_id + "." + imp.id if imp.id else file_id

        if not new_id:
            self.errors.set_file(file.path, file.name, self.options)
            self.error(
                imp.line, "No parent module -- cannot perform relative import", blocker=True
            )

        return new_id

    def all_imported_modules_in_file(self, file: MypyFile) -> list[tuple[int, str, int]]:
        """Find all reachable import statements in a file.

        Return list of tuples (priority, module id, import line number)
        for all modules imported in file; lower numbers == higher priority.

        Can generate blocking errors on bogus relative imports.
        """
        res: list[tuple[int, str, int]] = []
        for imp in file.imports:
            if not imp.is_unreachable:
                if isinstance(imp, Import):
                    pri = import_priority(imp, PRI_MED)
                    ancestor_pri = import_priority(imp, PRI_LOW)
                    for id, _ in imp.ids:
                        res.append((pri, id, imp.line))
                        ancestor_parts = id.split(".")[:-1]
                        ancestors = []
                        for part in ancestor_parts:
                            ancestors.append(part)
                            res.append((ancestor_pri, ".".join(ancestors), imp.line))
                elif isinstance(imp, ImportFrom):
                    cur_id = self.correct_rel_imp(file, imp)
                    all_are_submodules = True
                    # Also add any imported names that are submodules.
                    pri = import_priority(imp, PRI_MED)
                    for name, __ in imp.names:
                        sub_id = cur_id + "." + name
                        if self.is_module(sub_id):
                            res.append((pri, sub_id, imp.line))
                        else:
                            all_are_submodules = False
                    # Add cur_id as a dependency, even if all the
                    # imports are submodules. Processing import from will try
                    # to look through cur_id, so we should depend on it.
                    # As a workaround for some bugs in cycle handling (#4498),
                    # if all the imports are submodules, do the import at a lower
                    # priority.
                    pri = import_priority(imp, PRI_HIGH if not all_are_submodules else PRI_LOW)
                    res.append((pri, cur_id, imp.line))
                elif isinstance(imp, ImportAll):
                    pri = import_priority(imp, PRI_HIGH)
                    res.append((pri, self.correct_rel_imp(file, imp), imp.line))

        # Sort such that module (e.g. foo.bar.baz) comes before its ancestors (e.g. foo
        # and foo.bar) so that, if FindModuleCache finds the target module in a
        # package marked with py.typed underneath a namespace package installed in
        # site-packages, (gasp), that cache's knowledge of the ancestors
        # (aka FindModuleCache.ns_ancestors) can be primed when it is asked to find
        # the parent.
        res.sort(key=lambda x: -x[1].count("."))
        return res

    def is_module(self, id: str) -> bool:
        """Does the given fullname refer to a module?

        Note: this does not always verify that the module exists and relies on
        previously executed logic in find_module_and_diagnose().
        """
        if id in self.modules:
            # Micro-optimization, if we already found it, it is definitely a module.
            return True
        if id in self.source_set.source_modules:
            # Special case: if a module is passed on command line, we accept it even
            # if we would not resolve it using regular mechanisms. This makes behavior
            # consistent in cases like `mypy foo-stubs`, where stubs are not installed.
            return True
        return find_module_simple(id, self) is not None

    def parse_file(
        self,
        id: str,
        path: str,
        source: str,
        options: Options,
        raw_data: FileRawData | None = None,
    ) -> MypyFile:
        """Parse the source of a file with the given name.

        Raise CompileError if there is a parse error.
        """
        file_exists = self.fscache.exists(path, real_only=True)
        t0 = time.time()
        if raw_data:
            # If possible, deserialize from known binary data instead of parsing from scratch.
            tree = load_from_raw(path, id, raw_data, self.errors, options)
        else:
            tree = parse(source, path, id, self.errors, options=options, file_exists=file_exists)
        tree._fullname = id
        if self.stats_enabled:
            with self.stats_lock:
                self.add_stats(
                    files_parsed=1,
                    modules_parsed=int(not tree.is_stub),
                    stubs_parsed=int(tree.is_stub),
                    parse_time=time.time() - t0,
                )
        return tree

    def load_fine_grained_deps(self, id: str) -> dict[str, set[str]]:
        t0 = time.time()
        if id in self.fg_deps_meta:
            # TODO: Assert deps file wasn't changed.
            deps = json_loads(self.metastore.read(self.fg_deps_meta[id]["path"]))
        else:
            deps = {}
        val = {k: set(v) for k, v in deps.items()}
        self.add_stats(load_fg_deps_time=time.time() - t0)
        return val

    def report_file(
        self, file: MypyFile, type_map: dict[Expression, Type], options: Options
    ) -> None:
        if self.reports is not None and self.source_set.is_source(file):
            self.reports.file(file, self.modules, type_map, options)

    def commit(self) -> None:
        t0 = time.time()
        self.metastore.commit()
        self.add_stats(cache_commit_time=time.time() - t0)

    def commit_module(self, meta_file: str) -> None:
        """Commit cache writes for a single module (identified by its meta file path)."""
        self.metastore.commit_path(meta_file)

    def verbosity(self) -> int:
        return self.options.verbosity

    def log(self, *message: str) -> None:
        if self.verbosity() >= 1:
            if message:
                print("LOG: ", *message, file=self.stderr)
            else:
                print(file=self.stderr)
            self.stderr.flush()

    def log_fine_grained(self, *message: str) -> None:
        if self.verbosity() >= 1:
            self.log("fine-grained:", *message)
        elif mypy.build.DEBUG_FINE_GRAINED:
            # Output log in a simplified format that is quick to browse.
            if message:
                print(*message, file=self.stderr)
            else:
                print(file=self.stderr)
            self.stderr.flush()

    def trace(self, *message: str) -> None:
        if self.verbosity() >= 2:
            print("TRACE:", *message, file=self.stderr)
            self.stderr.flush()

    def add_stats(self, **kwds: Any) -> None:
        for key, value in kwds.items():
            if key in self.stats:
                self.stats[key] += value
            else:
                self.stats[key] = value

    def stats_summary(self) -> Mapping[str, object]:
        return self.stats

    def broadcast(self, message: bytes) -> None:
        """Broadcast same message to all workers in parallel."""
        t0 = time.time()
        threads = []
        for worker in self.workers:
            thread = Thread(target=worker.conn.write_bytes, args=(message,))
            thread.start()
            threads.append(thread)
        for thread in threads:
            thread.join()
        self.add_stats(broadcast_time=time.time() - t0)

    def wait_ack(self) -> None:
        """Wait for an ack from all workers."""
        for idx in range(len(self.workers)):
            buf = self.receive_worker_message(idx)
            assert read_tag(buf) == ACK_MESSAGE

    def receive_worker_message(self, idx: int) -> ReadBuffer:
        """Receive a single message from a worker, with crash diagnostics."""
        try:
            return receive(self.workers[idx].conn)
        except OSError as exc:
            try:
                # Give worker process a chance to actually terminate before reporting.
                exit_code = self.workers[idx].proc.wait(timeout=WORKER_SHUTDOWN_TIMEOUT)
            except TimeoutError:
                exit_code = None
            exit_status = f"exit code {exit_code}" if exit_code is not None else "still running"
            raise OSError(
                f"Worker {idx} disconnected before sending data ({exit_status})"
            ) from exc

    def submit(self, graph: Graph, sccs: list[SCC]) -> None:
        """Submit a stale SCC for processing in current process or parallel workers."""
        if self.workers:
            self.submit_to_workers(graph, sccs)
        else:
            self.scc_queue.extend([(0, 0, scc) for scc in sccs])

    def get_scc_batch(self, max_size_in_batch: int) -> list[SCC]:
        """Get a batch of SCCs from queue to submit to a worker.

        We batch SCCs to avoid communication overhead, but to avoid
        long poles, we limit fraction of work per worker.
        """
        batch: list[SCC] = []
        size_in_batch = 0
        while self.scc_queue and (
            # Three notes keep in mind here:
            #   * Heap key is *negative* size (so that larger SCCs appear first).
            #   * Each batch must have at least one item.
            #   * Adding another SCC to batch should not exceed maximum allowed size.
            size_in_batch - self.scc_queue[0][0] <= max_size_in_batch
            or not batch
        ):
            size_key, _, scc = heappop(self.scc_queue)
            size_in_batch -= size_key
            self.size_in_queue += size_key
            batch.append(scc)
        return batch

    def max_batch_size(self) -> int:
        batch_frac = 1 / len(self.workers)
        if sys.platform == "linux":
            # Linux is good with socket roundtrip latency, so we can use
            # more fine-grained batches.
            batch_frac /= 2
        return int(self.size_in_queue * batch_frac)

    def submit_to_workers(self, graph: Graph, sccs: list[SCC] | None = None) -> None:
        if sccs is not None:
            for scc in sccs:
                heappush(self.scc_queue, (-scc.size_hint, self.queue_order, scc))
                self.size_in_queue += scc.size_hint
                self.queue_order += 1
        max_size_in_batch = self.max_batch_size()
        while self.scc_queue and self.free_workers:
            idx = self.free_workers.pop()
            scc_batch = self.get_scc_batch(max_size_in_batch)
            import_errors = {
                mod_id: self.errors.recorded[path]
                for scc in scc_batch
                for mod_id in scc.mod_ids
                if (path := graph[mod_id].xpath) in self.errors.recorded
            }
            t0 = time.time()
            send(
                self.workers[idx].conn,
                SccRequestMessage(
                    scc_ids=[scc.id for scc in scc_batch],
                    import_errors=import_errors,
                    mod_data={
                        mod_id: (
                            # Although workers don't really need to know about details
                            # of dependencies, they will write cache, so we need to pass
                            # suppressed_deps_opts() as part of module data.
                            graph[mod_id].suppressed_deps_opts(),
                            tree.raw_data if (tree := graph[mod_id].tree) else None,
                        )
                        for scc in scc_batch
                        for mod_id in scc.mod_ids
                    },
                ),
            )
            self.add_stats(scc_requests_sent=1, scc_send_time=time.time() - t0)

    def wait_for_done(self, graph: Graph) -> tuple[list[SCC], bool, dict[str, ModuleResult]]:
        """Wait for a stale SCC processing to finish.

        Return a tuple three items:
          * processed SCCs
          * whether we have more in the queue
          * new interface hash or list of errors for each module
        The last item is only used for parallel processing.
        """
        if self.workers:
            return self.wait_for_done_workers(graph)
        if not self.scc_queue:
            return [], False, {}
        _, _, next_scc = self.scc_queue.pop(0)
        process_stale_scc(graph, next_scc, self)
        return [next_scc], bool(self.scc_queue), {}

    def wait_for_done_workers(
        self, graph: Graph
    ) -> tuple[list[SCC], bool, dict[str, ModuleResult]]:
        if not self.scc_queue and len(self.free_workers) == len(self.workers):
            return [], False, {}

        done_sccs = []
        results = {}
        t0 = time.time()
        ready = ready_to_read([w.conn for w in self.workers], WORKER_DONE_TIMEOUT)
        t1 = time.time()
        for idx in ready:
            buf = self.receive_worker_message(idx)
            assert read_tag(buf) == SCC_RESPONSE_MESSAGE
            data = SccResponseMessage.read(buf)
            if not data.is_interface:
                # Mark worker as free after it finished checking implementation.
                self.free_workers.add(idx)
            if data.blocker is not None:
                raise data.blocker
            assert data.result is not None
            results.update(data.result)
            if data.is_interface:
                done_sccs.extend([self.scc_by_id[scc_id] for scc_id in data.scc_ids])
        self.add_stats(scc_wait_time=t1 - t0, scc_receive_time=time.time() - t1)
        self.submit_to_workers(graph)  # advance after some workers are free.
        return (
            # Note that "done" means interface-ready in this context. This is what
            # the caller should expect.
            done_sccs,
            bool(self.scc_queue) or len(self.free_workers) < len(self.workers),
            results,
        )

    def is_transitive_scc_dep(self, from_scc_id: int, to_scc_id: int) -> bool:
        """Check if one SCC is a (transitive) dependency of another."""
        edge = (from_scc_id, to_scc_id)
        if (cached := self.transitive_deps_cache.get(edge)) is not None:
            return cached
        todo = self.scc_by_id[from_scc_id].deps
        seen = set()
        while todo:
            more = set()
            # Breadth-first search seems to be better here, because all
            # "lower-level" SCCs are processed and some may be cached.
            for dep in todo:
                seen.add(dep)
                if dep == to_scc_id:
                    self.transitive_deps_cache[edge] = True
                    return True
                if cached := self.transitive_deps_cache.get((dep, to_scc_id)):
                    self.transitive_deps_cache[edge] = True
                    return True
                elif cached is None:
                    more |= self.scc_by_id[dep].deps
            todo = more
        self.transitive_deps_cache[edge] = False
        for dep in seen:
            # We negative-cache all intermediate lookups, thus
            # trading time for space.
            self.transitive_deps_cache[(dep, to_scc_id)] = False
        return False

    def error(
        self,
        line: int | None,
        msg: str,
        code: ErrorCode | None = None,
        *,
        blocker: bool = False,
        only_once: bool = False,
    ) -> None:
        if line is None:
            line = column = -1
        else:
            column = 0
        self.errors.report(line, column, msg, code, blocker=blocker, only_once=only_once)

    def note(
        self, line: int | None, msg: str, code: ErrorCode | None = None, *, only_once: bool = False
    ) -> None:
        if line is None:
            line = column = -1
        else:
            column = 0
        self.errors.report(line, column, msg, code, severity="note", only_once=only_once)

    def note_multiline(
        self, line: int | None, msg: str, code: ErrorCode | None = None, *, only_once: bool = False
    ) -> None:
        for msg_line in dedent(msg.lstrip("\n")).splitlines():
            self.note(line, msg_line, code, only_once=only_once)


def deps_to_json(x: dict[str, set[str]]) -> bytes:
    return json_dumps({k: list(v) for k, v in x.items()})


# File for storing metadata about all the fine-grained dependency caches
DEPS_META_FILE: Final = "@deps.meta.json"
# File for storing fine-grained dependencies that didn't a parent in the build
DEPS_ROOT_FILE: Final = "@root.deps.json"

# The name of the fake module used to store fine-grained dependencies that
# have no other place to go.
FAKE_ROOT_MODULE: Final = "@root"


def write_deps_cache(
    rdeps: dict[str, dict[str, set[str]]], manager: BuildManager, graph: Graph
) -> None:
    """Write cache files for fine-grained dependencies.

    Serialize fine-grained dependencies map for fine-grained mode.

    Dependencies on some module 'm' is stored in the dependency cache
    file m.deps.json.  This entails some spooky action at a distance:
    if module 'n' depends on 'm', that produces entries in m.deps.json.
    When there is a dependency on a module that does not exist in the
    build, it is stored with its first existing parent module. If no
    such module exists, it is stored with the fake module FAKE_ROOT_MODULE.

    This means that the validity of the fine-grained dependency caches
    are a global property, so we store validity checking information for
    fine-grained dependencies in a global cache file:
     * We take a snapshot of current sources to later check consistency
       between the fine-grained dependency cache and module cache metadata
     * We store the mtime of all the dependency files to verify they
       haven't changed
    """
    metastore = manager.metastore

    error = False

    fg_deps_meta = manager.fg_deps_meta.copy()

    for id in rdeps:
        if id != FAKE_ROOT_MODULE:
            _, _, deps_json = get_cache_names(id, graph[id].xpath, manager.options)
        else:
            deps_json = DEPS_ROOT_FILE
        assert deps_json
        manager.log("Writing deps cache", deps_json)
        if not manager.metastore.write(deps_json, deps_to_json(rdeps[id])):
            manager.log(f"Error writing fine-grained deps JSON file {deps_json}")
            error = True
        else:
            fg_deps_meta[id] = {"path": deps_json, "mtime": manager.getmtime(deps_json)}

    meta_snapshot: dict[str, str] = {}
    for id, st in graph.items():
        # If we didn't parse a file (so it doesn't have a
        # source_hash), then it must be a module with a fresh cache,
        # so use the hash from that.
        if st.source_hash:
            hash = st.source_hash
        else:
            if st.meta:
                hash = st.meta.hash
            else:
                hash = ""
        meta_snapshot[id] = hash

    meta = {"snapshot": meta_snapshot, "deps_meta": fg_deps_meta}

    if not metastore.write(DEPS_META_FILE, json_dumps(meta)):
        manager.log(f"Error writing fine-grained deps meta JSON file {DEPS_META_FILE}")
        error = True

    if error:
        manager.errors.set_file(_cache_dir_prefix(manager.options), None, manager.options)
        manager.error(None, "Error writing fine-grained dependencies cache", blocker=True)


def invert_deps(deps: dict[str, set[str]], graph: Graph) -> dict[str, dict[str, set[str]]]:
    """Splits fine-grained dependencies based on the module of the trigger.

    Returns a dictionary from module ids to all dependencies on that
    module. Dependencies not associated with a module in the build will be
    associated with the nearest parent module that is in the build, or the
    fake module FAKE_ROOT_MODULE if none are.
    """
    # Lazy import to speed up startup
    from mypy.server.target import trigger_to_target

    # Prepopulate the map for all the modules that have been processed,
    # so that we always generate files for processed modules (even if
    # there aren't any dependencies to them.)
    rdeps: dict[str, dict[str, set[str]]] = {id: {} for id, st in graph.items() if st.tree}
    for trigger, targets in deps.items():
        module = module_prefix(graph, trigger_to_target(trigger))
        if not module or not graph[module].tree:
            module = FAKE_ROOT_MODULE

        mod_rdeps = rdeps.setdefault(module, {})
        mod_rdeps.setdefault(trigger, set()).update(targets)

    return rdeps


def generate_deps_for_cache(manager: BuildManager, graph: Graph) -> dict[str, dict[str, set[str]]]:
    """Generate fine-grained dependencies into a form suitable for serializing.

    This does a couple things:
    1. Splits fine-grained deps based on the module of the trigger
    2. For each module we generated fine-grained deps for, load any previous
       deps and merge them in.

    Returns a dictionary from module ids to all dependencies on that
    module. Dependencies not associated with a module in the build will be
    associated with the nearest parent module that is in the build, or the
    fake module FAKE_ROOT_MODULE if none are.
    """
    from mypy.server.deps import merge_dependencies  # Lazy import to speed up startup

    # Split the dependencies out into based on the module that is depended on.
    rdeps = invert_deps(manager.fg_deps, graph)

    # We can't just clobber existing dependency information, so we
    # load the deps for every module we've generated new dependencies
    # to and merge the new deps into them.
    for module, mdeps in rdeps.items():
        old_deps = manager.load_fine_grained_deps(module)
        merge_dependencies(old_deps, mdeps)

    return rdeps


PLUGIN_SNAPSHOT_FILE: Final = "@plugins_snapshot.json"


def write_plugins_snapshot(manager: BuildManager) -> None:
    """Write snapshot of versions and hashes of currently active plugins."""
    snapshot = json_dumps(manager.plugins_snapshot)
    if (
        not manager.metastore.write(PLUGIN_SNAPSHOT_FILE, snapshot)
        and manager.options.cache_dir != os.devnull
    ):
        manager.errors.set_file(_cache_dir_prefix(manager.options), None, manager.options)
        manager.error(None, "Error writing plugins snapshot", blocker=True)


def read_plugins_snapshot(manager: BuildManager) -> dict[str, str] | None:
    """Read cached snapshot of versions and hashes of plugins from previous run."""
    snapshot = _load_json_file(
        PLUGIN_SNAPSHOT_FILE,
        manager,
        log_success="Plugins snapshot (cached) ",
        log_error="Could not load plugins snapshot: ",
    )
    if snapshot is None:
        return None
    if not isinstance(snapshot, dict):
        manager.log(f"Could not load plugins snapshot: cache is not a dict: {type(snapshot)}")  # type: ignore[unreachable]
        return None
    return snapshot


def read_quickstart_file(
    options: Options, stdout: TextIO
) -> dict[str, tuple[float, int, str]] | None:
    quickstart: dict[str, tuple[float, int, str]] | None = None
    if options.quickstart_file:
        # This is very "best effort". If the file is missing or malformed,
        # just ignore it.
        raw_quickstart: dict[str, Any] = {}
        try:
            with open(options.quickstart_file, "rb") as f:
                raw_quickstart = json_loads(f.read())

            quickstart = {}
            for file, (x, y, z) in raw_quickstart.items():
                quickstart[file] = (x, y, z)
        except Exception as e:
            print(f"Warning: Failed to load quickstart file: {str(e)}\n", file=stdout)
    return quickstart


def read_deps_cache(manager: BuildManager, graph: Graph) -> dict[str, FgDepMeta] | None:
    """Read and validate the fine-grained dependencies cache.

    See the write_deps_cache documentation for more information on
    the details of the cache.

    Returns None if the cache was invalid in some way.
    """
    deps_meta = _load_json_file(
        DEPS_META_FILE,
        manager,
        log_success="Deps meta ",
        log_error="Could not load fine-grained dependency metadata: ",
    )
    if deps_meta is None:
        return None
    meta_snapshot = deps_meta["snapshot"]
    # Take a snapshot of the source hashes from all the metas we found.
    # (Including the ones we rejected because they were out of date.)
    # We use this to verify that they match up with the proto_deps.
    current_meta_snapshot = {
        id: st.meta_source_hash for id, st in graph.items() if st.meta_source_hash is not None
    }

    common = set(meta_snapshot.keys()) & set(current_meta_snapshot.keys())
    if any(meta_snapshot[id] != current_meta_snapshot[id] for id in common):
        # TODO: invalidate also if options changed (like --strict-optional)?
        manager.log("Fine-grained dependencies cache inconsistent, ignoring")
        return None

    module_deps_metas = deps_meta["deps_meta"]
    assert isinstance(module_deps_metas, dict)
    if not manager.options.skip_cache_mtime_checks:
        for meta in module_deps_metas.values():
            try:
                matched = manager.getmtime(meta["path"]) == meta["mtime"]
            except FileNotFoundError:
                matched = False
            if not matched:
                manager.log(f"Invalid or missing fine-grained deps cache: {meta['path']}")
                return None

    return module_deps_metas


def _load_ff_file(
    file: str, manager: BuildManager, log_error_fmt: str, id: str | None
) -> bytes | None:
    if manager.stats_enabled:
        t0 = time.time()
    try:
        data = manager.metastore.read(file)
    except OSError:
        if manager.logging_enabled:
            if id:
                message = log_error_fmt.format(id) + file
            else:
                message = log_error_fmt + file
            manager.log(message)
        return None
    if manager.stats_enabled:
        manager.add_stats(metastore_read_time=time.time() - t0)
    return data


def _load_json_file(
    file: str, manager: BuildManager, log_success: str, log_error: str
) -> dict[str, Any] | None:
    """A simple helper to read a JSON file with logging."""
    t0 = time.time()
    try:
        data = manager.metastore.read(file)
    except OSError:
        manager.log(log_error + file)
        return None
    manager.add_stats(metastore_read_time=time.time() - t0)
    # Only bother to compute the log message if we are logging it, since it could be big
    if manager.verbosity() >= 2:
        manager.trace(log_success + data.rstrip().decode())
    try:
        t1 = time.time()
        result = json_loads(data)
        manager.add_stats(data_file_load_time=time.time() - t1)
    except json.JSONDecodeError:
        manager.errors.set_file(file, None, manager.options)
        manager.error(
            None,
            "Error reading JSON file;"
            " you likely have a bad cache.\n"
            "Try removing the {cache_dir} directory"
            " and run mypy again.".format(cache_dir=manager.options.cache_dir),
            blocker=True,
        )
        return None
    else:
        assert isinstance(result, dict)
        return result


def _cache_dir_prefix(options: Options) -> str:
    """Get current cache directory (or file if id is given)."""
    if options.bazel:
        # This is needed so the cache map works.
        return os.curdir
    cache_dir = options.cache_dir
    pyversion = options.python_version
    base = os_path_join(cache_dir, "%d.%d" % pyversion)
    return base


def add_catch_all_gitignore(target_dir: str) -> None:
    """Add catch-all .gitignore to an existing directory.

    No-op if the .gitignore already exists.
    """
    gitignore = os_path_join(target_dir, ".gitignore")
    try:
        with open(gitignore, "x") as f:
            print("# Automatically created by mypy", file=f)
            print("*", file=f)
    except FileExistsError:
        pass


def exclude_from_backups(target_dir: str) -> None:
    """Exclude the directory from various archives and backups supporting CACHEDIR.TAG.

    If the CACHEDIR.TAG file exists the function is a no-op.
    """
    cachedir_tag = os_path_join(target_dir, "CACHEDIR.TAG")
    try:
        with open(cachedir_tag, "x") as f:
            f.write("""Signature: 8a477f597d28d172789f06886806bc55
# This file is a cache directory tag automatically created by mypy.
# For information about cache directory tags see https://bford.info/cachedir/
""")
    except FileExistsError:
        pass


def create_metastore(options: Options, parallel_worker: bool) -> MetadataStore:
    """Create the appropriate metadata store."""
    if options.sqlite_cache:
        mds: MetadataStore = SqliteMetadataStore(
            _cache_dir_prefix(options),
            set_journal_mode=not parallel_worker,
            num_shards=options.sqlite_num_shards,
        )
    else:
        mds = FilesystemMetadataStore(_cache_dir_prefix(options))
    return mds


def get_meta_ex_name(meta_name: str) -> str:
    # Convert e.g. foo.bar.meta.ff to foo.bar.meta_ex.ff
    parts = meta_name.rsplit(".", maxsplit=2)
    parts[1] = "meta_ex"
    return ".".join(parts)


def get_cache_names(id: str, path: str, options: Options) -> tuple[str, str, str | None]:
    """Return the file names for the cache files.

    Args:
      id: module ID
      path: module path
      options: build options

    Returns:
      A tuple with the file names to be used for the meta file, the
      data file, and the fine-grained deps JSON, respectively.
    """
    if options.cache_map:
        pair = options.cache_map.get(normpath(path, options))
    else:
        pair = None
    if pair is not None:
        # The cache map paths were specified relative to the base directory,
        # but the filesystem metastore APIs operates relative to the cache
        # prefix directory.
        # Solve this by rewriting the paths as relative to the root dir.
        # This only makes sense when using the filesystem backed cache.
        root = _cache_dir_prefix(options)
        return os.path.relpath(pair[0], root), os.path.relpath(pair[1], root), None
    prefix = os.path.join(*id.split("."))
    is_package = os.path.basename(path).startswith("__init__.py")
    if is_package:
        prefix = os_path_join(prefix, "__init__")

    deps_json = None
    if options.cache_fine_grained:
        deps_json = prefix + ".deps.json"
    if options.fixed_format_cache:
        data_suffix = ".data.ff"
        meta_suffix = ".meta.ff"
    else:
        data_suffix = ".data.json"
        meta_suffix = ".meta.json"
    return prefix + meta_suffix, prefix + data_suffix, deps_json


def options_snapshot(id: str, manager: BuildManager) -> dict[str, object]:
    """Make compact snapshot of options for a module.

    Separately store only the options we may compare individually, and take a hash
    of everything else. If --debug-cache is specified, fall back to full snapshot.
    """
    platform_opt, values = manager.options.clone_for_module(id).select_options_affecting_cache()
    if manager.options.debug_cache:
        # Build full options snapshot for debugging purposes.
        result: dict[str, object] = {"platform": platform_opt}
        for key, val in zip(OPTIONS_AFFECTING_CACHE_NO_PLATFORM, values):
            result[key] = val
        return result
    # Process most options quickly, since this is performance critical.
    buf = WriteBuffer()
    write_json_value(buf, cast(JsonValue, values))
    return {"platform": platform_opt, "other_options": hash_digest(buf.getvalue())}


def find_cache_meta(
    id: str, path: str, manager: BuildManager, skip_validation: bool = False
) -> tuple[CacheMeta, CacheMetaEx] | None:
    """Find cache data for a module.

    Args:
      id: module ID
      path: module path
      manager: the build manager (for pyversion, log/trace, and build options)
      skip_validation: if True skip any validation steps (used for parallel checking)

    Returns:
      A CacheMeta/CacheMetaEx instance pair if the cache data was found and appears
      valid; otherwise None.
    """
    # TODO: May need to take more build options into account
    meta_file, data_file, _ = get_cache_names(id, path, manager.options)
    if manager.tracing_enabled:
        manager.trace(f"Looking for {id} at {meta_file}")
    if manager.stats_enabled:
        t0 = time.time()
    if manager.options.fixed_format_cache:
        meta = _load_ff_file(
            meta_file, manager, log_error_fmt="Could not load cache for {}: ", id=id
        )
    else:
        meta = _load_json_file(
            meta_file,
            manager,
            log_success=f"Meta {id} ",
            log_error=f"Could not load cache for {id}: ",
        )
    if meta is None:
        return None
    if manager.stats_enabled:
        t1 = time.time()
    if isinstance(meta, bytes):
        # If either low-level buffer format or high-level cache layout changed, we
        # cannot use the cache files, even with --skip-version-check.
        # TODO: switch to something like librt.internal.read_byte() if this is slow.
        if meta[0] != cache_version() or meta[1] != CACHE_VERSION:
            manager.log(f"Metadata abandoned for {id}: incompatible cache format")
            return None
        data_io = ReadBuffer(meta[2:])
        m = CacheMeta.read(data_io, data_file)
    else:
        m = CacheMeta.deserialize(meta, data_file)
    if m is None:
        manager.log(f"Metadata abandoned for {id}: cannot deserialize data")
        return None
    if manager.stats_enabled:
        t2 = time.time()
        manager.add_stats(
            load_meta_time=t2 - t0, load_meta_load_time=t1 - t0, load_meta_from_dict_time=t2 - t1
        )
    if skip_validation:
        # If the caller requested no validation, skip the implementation part of the meta
        # as well, as a performance optimization. Note: this may return an incomplete meta,
        # only use if you know what you are doing.
        assert manager.parallel_worker
        return m, CacheMetaEx([], [], [], [])

    # Ignore cache if generated by an older mypy version.
    if m.version_id != manager.version_id and not manager.options.skip_version_check:
        manager.log(f"Metadata abandoned for {id}: different mypy version")
        return None

    total_deps = len(m.dependencies) + len(m.suppressed)
    if len(m.dep_prios) != total_deps or len(m.dep_lines) != total_deps:
        manager.log(f"Metadata abandoned for {id}: broken dependencies")
        return None

    # Ignore cache if (relevant) options aren't the same.
    # Note that it's fine to mutilate cached_options since it's only used here.
    cached_options = m.options
    current_options = options_snapshot(id, manager)
    if manager.options.skip_version_check:
        # When we're lax about version we're also lax about platform.
        cached_options["platform"] = current_options["platform"]
    if "debug_cache" in cached_options:
        # Older versions included debug_cache, but it's silly to compare it.
        del cached_options["debug_cache"]
    if cached_options != current_options:
        manager.log(f"Metadata abandoned for {id}: options differ")
        if manager.options.verbosity >= 2:
            for key in sorted(set(cached_options) | set(current_options)):
                if cached_options.get(key) != current_options.get(key):
                    manager.trace(
                        "    {}: {} != {}".format(
                            key, cached_options.get(key), current_options.get(key)
                        )
                    )
        return None
    if manager.old_plugins_snapshot and manager.plugins_snapshot:
        # Check if plugins are still the same.
        if manager.plugins_snapshot != manager.old_plugins_snapshot:
            manager.log(f"Metadata abandoned for {id}: plugins differ")
            return None
    plugin_data = manager.plugin.report_config_data(ReportConfigContext(id, path, is_check=True))
    if not manager.options.fixed_format_cache:
        # So that plugins can return data with tuples in it without
        # things silently always invalidating modules, we round-trip
        # the config data. This isn't beautiful.
        plugin_data = json_loads(json_dumps(plugin_data))
    if m.plugin_data != plugin_data:
        manager.log(f"Metadata abandoned for {id}: plugin configuration differs")
        return None

    meta_ex_file = get_meta_ex_name(meta_file)
    if manager.options.fixed_format_cache:
        meta_ex = _load_ff_file(
            meta_ex_file, manager, log_error_fmt="Could not load meta_ex for {}: ", id=id
        )
    else:
        meta_ex = _load_json_file(
            meta_ex_file,
            manager,
            log_success=f"Meta_ex {id} ",
            log_error=f"Could not load meta_ex for {id}: ",
        )
    if meta_ex is None:
        return None
    if isinstance(meta_ex, bytes):
        data_io = ReadBuffer(meta_ex)
        me = CacheMetaEx.read(data_io)
    else:
        me = CacheMetaEx.deserialize(meta_ex)
    if me is None:
        return None
    manager.add_stats(fresh_metas=1)
    return m, me


def validate_meta(
    meta: CacheMeta | None, id: str, path: str | None, ignore_all: bool, manager: BuildManager
) -> CacheMeta | None:
    """Checks whether the cached AST of this module can be used.

    Returns:
      None, if the cached AST is unusable.
      Original meta, if mtime/size matched.
      Meta with mtime updated to match source file, if hash/size matched but mtime/path didn't.
    """
    # This requires two steps. The first one is obvious: we check that the module source file
    # contents is the same as it was when the cache data file was created. The second one is not
    # too obvious: we check that the cache data file mtime has not changed; it is needed because
    # we use cache data file mtime to propagate information about changes in the dependencies.

    if meta is None:
        manager.log(f"Metadata not found for {id}")
        return None

    if meta.ignore_all and not ignore_all:
        manager.log(f"Metadata abandoned for {id}: errors were previously ignored")
        return None

    if manager.stats_enabled:
        t0 = time.time()
    bazel = manager.options.bazel
    assert path is not None, "Internal error: meta was provided without a path"
    if not manager.options.skip_cache_mtime_checks:
        # Check data_file; assume if its mtime matches it's good.
        try:
            data_mtime = manager.getmtime(meta.data_file)
        except OSError:
            manager.log(f"Metadata abandoned for {id}: failed to stat data_file")
            return None
        if data_mtime != meta.data_mtime:
            manager.log(f"Metadata abandoned for {id}: data cache is modified")
            return None

    if bazel:
        # Normalize path under bazel to make sure it isn't absolute
        path = normpath(path, manager.options)

    st = manager.get_stat(path)
    if st is None:
        return None
    if not stat.S_ISDIR(st.st_mode) and not stat.S_ISREG(st.st_mode):
        manager.log(f"Metadata abandoned for {id}: file or directory {path} does not exist")
        return None

    if manager.stats_enabled:
        manager.add_stats(validate_stat_time=time.time() - t0)

    # When we are using a fine-grained cache, we want our initial
    # build() to load all of the cache information and then do a
    # fine-grained incremental update to catch anything that has
    # changed since the cache was generated. We *don't* want to do a
    # coarse-grained incremental rebuild, so we accept the cache
    # metadata even if it doesn't match the source file.
    #
    # We still *do* the mtime/hash checks, however, to enable
    # fine-grained mode to take advantage of the mtime-updating
    # optimization when mtimes differ but hashes match.  There is
    # essentially no extra time cost to computing the hash here, since
    # it will be cached and will be needed for finding changed files
    # later anyways.
    fine_grained_cache = manager.use_fine_grained_cache()

    size = st.st_size
    # Bazel ensures the cache is valid.
    if size != meta.size and not bazel and not fine_grained_cache:
        manager.log(f"Metadata abandoned for {id}: file {path} has different size")
        return None

    # Bazel ensures the cache is valid.
    mtime = 0 if bazel else int(st.st_mtime)
    if not bazel and (mtime != meta.mtime or path != meta.path):
        if manager.quickstart_state and path in manager.quickstart_state:
            # If the mtime and the size of the file recorded in the quickstart dump matches
            # what we see on disk, we know (assume) that the hash matches the quickstart
            # data as well. If that hash matches the hash in the metadata, then we know
            # the file is up to date even though the mtime is wrong, without needing to hash it.
            qmtime, qsize, qhash = manager.quickstart_state[path]
            if int(qmtime) == mtime and qsize == size and qhash == meta.hash:
                manager.log(f"Metadata fresh (by quickstart) for {id}: file {path}")
                meta.mtime = mtime
                meta.path = path
                return meta

        t0 = time.time()
        try:
            # dir means it is a namespace package
            if stat.S_ISDIR(st.st_mode):
                source_hash = ""
            else:
                source_hash = manager.fscache.hash_digest(path)
        except (OSError, UnicodeDecodeError, DecodeError):
            return None
        manager.add_stats(validate_hash_time=time.time() - t0)
        if source_hash != meta.hash:
            if fine_grained_cache:
                manager.log(f"Using stale metadata for {id}: file {path}")
                return meta
            else:
                manager.log(f"Metadata abandoned for {id}: file {path} has different hash")
                return None
        else:
            if manager.stats_enabled:
                t0 = time.time()
            # Optimization: update mtime and path (otherwise, this mismatch will reappear).
            meta.mtime = mtime
            meta.path = path
            meta.size = size
            meta.options = options_snapshot(id, manager)
            meta_file, _, _ = get_cache_names(id, path, manager.options)
            if manager.logging_enabled:
                manager.log(
                    "Updating mtime for {}: file {}, meta {}, mtime {}".format(
                        id, path, meta_file, meta.mtime
                    )
                )
            write_cache_meta(meta, manager, meta_file)
            if manager.stats_enabled:
                t1 = time.time()
                manager.add_stats(
                    validate_update_time=time.time() - t1, validate_munging_time=t1 - t0
                )
            return meta

    # It's a match on (id, path, size, hash, mtime).
    if manager.logging_enabled:
        manager.log(f"Metadata fresh for {id}: file {path}")
    return meta


def compute_hash(text: str) -> str:
    # We use a crypto hash instead of the builtin hash(...) function
    # because the output of hash(...)  can differ between runs due to
    # hash randomization (enabled by default in Python 3.3).  See the
    # note in
    # https://docs.python.org/3/reference/datamodel.html#object.__hash__.
    return hash_digest(text.encode("utf-8"))


def write_cache(
    id: str,
    path: str,
    tree: MypyFile,
    dependencies: list[str],
    suppressed: list[str],
    suppressed_deps_opts: bytes,
    imports_ignored: dict[int, list[str]],
    dep_prios: list[int],
    dep_lines: list[int],
    old_interface_hash: bytes,
    trans_dep_hash: bytes,
    source_hash: str,
    ignore_all: bool,
    manager: BuildManager,
) -> tuple[bytes, tuple[CacheMeta, str] | None]:
    """Write cache files for a module.

    Note that this mypy's behavior is still correct when any given
    write_cache() call is replaced with a no-op, so error handling
    code that bails without writing anything is okay.

    Args:
      id: module ID
      path: module path
      tree: the fully checked module data
      dependencies: module IDs on which this module depends
      suppressed: module IDs which were suppressed as dependencies
      dep_prios: priorities (parallel array to dependencies)
      dep_lines: import line locations (parallel array to dependencies)
      old_interface_hash: the hash from the previous version of the data cache file
      source_hash: the hash of the source code
      ignore_all: the ignore_all flag for this module
      manager: the build manager (for pyversion, log/trace)

    Returns:
      A tuple containing the interface hash and inner tuple with CacheMeta
      that should be written and path to cache file (inner tuple may be None,
      if the cache data could not be written).
    """
    metastore = manager.metastore
    # For Bazel we use relative paths and zero mtimes.
    bazel = manager.options.bazel

    # Obtain file paths.
    meta_file, data_file, _ = get_cache_names(id, path, manager.options)
    manager.log(f"Writing {id} {path} {meta_file} {data_file}")

    # Update tree.path so that in bazel mode it's made relative (since
    # sometimes paths leak out).
    if bazel:
        tree.path = path

    plugin_data = manager.plugin.report_config_data(ReportConfigContext(id, path, is_check=False))

    # Serialize data and analyze interface
    if manager.options.fixed_format_cache:
        data_io = WriteBuffer()
        tree.write(data_io)
        data_bytes = data_io.getvalue()
    else:
        data = tree.serialize()
        data_bytes = json_dumps(data, manager.options.debug_cache)
    interface_hash = hash_digest_bytes(data_bytes + json_dumps(plugin_data))

    # Obtain and set up metadata
    st = manager.get_stat(path)
    if st is None:
        manager.log(f"Cannot get stat for {path}")
        # Remove apparently-invalid cache files.
        # (This is purely an optimization.)
        for filename in [data_file, meta_file]:
            try:
                os.remove(filename)
            except OSError:
                pass
        # Still return the interface hash we computed.
        return interface_hash, None

    # Write data cache file, if applicable
    # Note that for Bazel we don't record the data file's mtime.
    if old_interface_hash == interface_hash:
        manager.trace(f"Interface for {id} is unchanged")
    else:
        manager.trace(f"Interface for {id} has changed")
        if not metastore.write(data_file, data_bytes):
            # Most likely the error is the replace() call
            # (see https://github.com/python/mypy/issues/3215).
            manager.log(f"Error writing cache data file {data_file}")
            # Let's continue without writing the meta file.  Analysis:
            # If the replace failed, we've changed nothing except left
            # behind an extraneous temporary file; if the replace
            # worked but the getmtime() call failed, the meta file
            # will be considered invalid on the next run because the
            # data_mtime field won't match the data file's mtime.
            # Both have the effect of slowing down the next run a
            # little bit due to an out-of-date cache file.
            return interface_hash, None

    try:
        data_mtime = manager.getmtime(data_file)
    except OSError:
        manager.log(f"Error in os.stat({data_file!r}), skipping cache write")
        return interface_hash, None

    mtime = 0 if bazel else int(st.st_mtime)
    size = st.st_size
    # Note that the options we store in the cache are the options as
    # specified by the command line/config file and *don't* reflect
    # updates made by inline config directives in the file. This is
    # important, or otherwise the options would never match when
    # verifying the cache.
    assert source_hash is not None
    meta = CacheMeta(
        id=id,
        path=path,
        mtime=mtime,
        size=size,
        hash=source_hash,
        dependencies=dependencies,
        data_mtime=data_mtime,
        data_file=data_file,
        suppressed=suppressed,
        imports_ignored=imports_ignored,
        options=options_snapshot(id, manager),
        suppressed_deps_opts=suppressed_deps_opts,
        dep_prios=dep_prios,
        dep_lines=dep_lines,
        interface_hash=interface_hash,
        trans_dep_hash=trans_dep_hash,
        version_id=manager.version_id,
        ignore_all=ignore_all,
        plugin_data=plugin_data,
        # This one will be filled by the caller.
        dep_hashes=[],
    )
    return interface_hash, (meta, meta_file)


def write_cache_meta(meta: CacheMeta, manager: BuildManager, meta_file: str) -> None:
    # Write meta cache file
    metastore = manager.metastore
    if manager.options.fixed_format_cache:
        data_io = WriteBuffer()
        meta.write(data_io)
        # Prefix with both low- and high-level cache format versions for future validation.
        # TODO: switch to something like librt.internal.write_byte() if this is slow.
        meta_bytes = bytes([cache_version(), CACHE_VERSION]) + data_io.getvalue()
    else:
        meta_dict = meta.serialize()
        meta_bytes = json_dumps(meta_dict, manager.options.debug_cache)
    if not metastore.write(meta_file, meta_bytes):
        # Most likely the error is the replace() call
        # (see https://github.com/python/mypy/issues/3215).
        # The next run will simply find the cache entry out of date.
        manager.log(f"Error writing cache meta file {meta_file}")


def write_cache_meta_ex(meta_file: str, meta_ex: CacheMetaEx, manager: BuildManager) -> None:
    # Write errors cache file
    meta_ex_file = get_meta_ex_name(meta_file)
    metastore = manager.metastore
    if manager.options.fixed_format_cache:
        data_io = WriteBuffer()
        meta_ex.write(data_io)
        meta_bytes = data_io.getvalue()
    else:
        # Some generic JSON helpers require top-level to be a dict.
        meta_bytes = json_dumps(meta_ex.serialize(), manager.options.debug_cache)
    if not metastore.write(meta_ex_file, meta_bytes):
        manager.log(f"Error writing meta_ex file {meta_ex_file}")


"""Dependency manager.

Design
======

Ideally
-------

A. Collapse cycles (each SCC -- strongly connected component --
   becomes one "supernode").

B. Topologically sort nodes based on dependencies.

C. Process from leaves towards roots.

Wrinkles
--------

a. Need to parse source modules to determine dependencies.

b. Processing order for modules within an SCC.

c. Must order mtimes of files to decide whether to re-process; depends
   on clock never resetting.

d. from P import M; checks filesystem whether module P.M exists in
   filesystem.

e. Race conditions, where somebody modifies a file while we're
   processing. Solved by using a FileSystemCache.


Steps
-----

1. For each explicitly given module find the source file location.

2. For each such module load and check the cache metadata, and decide
   whether it's valid.

3. Now recursively (or iteratively) find dependencies and add those to
   the graph:

   - for cached nodes use the list of dependencies from the cache
     metadata (this will be valid even if we later end up re-parsing
     the same source);

   - for uncached nodes parse the file and process all imports found,
     taking care of (a) above.

Step 3 should also address (d) above.

Once step 3 terminates we have the entire dependency graph, and for
each module we've either loaded the cache metadata or parsed the
source code.  (However, we may still need to parse those modules for
which we have cache metadata but that depend, directly or indirectly,
on at least one module for which the cache metadata is stale.)

Now we can execute steps A-C from the first section.  Finding SCCs for
step A shouldn't be hard; there's a recipe here:
https://code.activestate.com/recipes/578507/.  There's also a plethora
of topsort recipes, e.g. https://code.activestate.com/recipes/577413/.

For single nodes, processing is simple.  If the node was cached, we
deserialize the cache data and fix up cross-references.  Otherwise, we
do semantic analysis followed by type checking.  Once we (re-)processed
an SCC we check whether its interface (symbol table) is still fresh
(matches previous cached value). If it is not, we consider dependent SCCs
stale so that they need to be re-parsed as well.

Note on indirect dependencies: normally dependencies are determined from
imports, but since our interfaces are "opaque" (i.e. symbol tables can
contain cross-references as well as types identified by name), these are not
enough.  We *must* also add "indirect" dependencies from symbols and types to
their definitions.  For this purpose, we record all accessed symbols during
semantic analysis, and after we finished processing a module, we traverse its
type map, and for each type we find (transitively) on which named types it
depends.

Import cycles
-------------

Finally we have to decide how to handle (b), import cycles.  Here
we'll need a modified version of the original state machine
(build.py), but we only need to do this per SCC, and we won't have to
deal with changes to the list of nodes while we're processing it.

If all nodes in the SCC have valid cache metadata and all dependencies
outside the SCC are still valid, we can proceed as follows:

  1. Load cache data for all nodes in the SCC.

  2. Fix up cross-references for all nodes in the SCC.

Otherwise, the simplest (but potentially slow) way to proceed is to
invalidate all cache data in the SCC and re-parse all nodes in the SCC
from source.  We can do this as follows:

  1. Parse source for all nodes in the SCC.

  2. Semantic analysis for all nodes in the SCC.

  3. Type check all nodes in the SCC.

(If there are more passes the process is the same -- each pass should
be done for all nodes before starting the next pass for any nodes in
the SCC.)

We could process the nodes in the SCC in any order.  For sentimental
reasons, I've decided to process them in the reverse order in which we
encountered them when originally constructing the graph.  That's how
the old build.py deals with cycles, and at least this reproduces the
previous implementation more accurately.

Can we do better than re-parsing all nodes in the SCC when any of its
dependencies are out of date?  It's doubtful.  The optimization
mentioned at the end of the previous section would require re-parsing
and type-checking a node and then comparing its symbol table to the
cached data; but because the node is part of a cycle we can't
technically type-check it until the semantic analysis of all other
nodes in the cycle has completed.  (This is an important issue because
Dropbox has a very large cycle in production code.  But I'd like to
deal with it later.)

Additional wrinkles
-------------------

During implementation more wrinkles were found.

- When a submodule of a package (e.g. x.y) is encountered, the parent
  package (e.g. x) must also be loaded, but it is not strictly a
  dependency.  See State.add_ancestors() below.
"""


class SuppressionReason:
    NOT_FOUND: Final = 1
    SKIPPED: Final = 2


class ModuleNotFound(Exception):
    """Control flow exception to signal that a module was not found."""

    def __init__(self, reason: int = SuppressionReason.NOT_FOUND) -> None:
        self.reason = reason


@final
class State:
    """The state for a module.

    The source is only used for the -c command line option; in that
    case path is None.  Otherwise, source is None and path isn't.
    """

    manager: BuildManager
    order_counter: ClassVar[int] = 0
    order: int  # Order in which modules were encountered
    id: str  # Fully qualified module name
    path: str | None = None  # Path to module source
    abspath: str | None = None  # Absolute path to module source
    xpath: str  # Path or '<string>'
    source: str | None = None  # Module source code
    source_hash: str | None = None  # Hash calculated based on the source code
    meta_source_hash: str | None = None  # Hash of the source given in the meta, if any
    meta: CacheMeta | None = None
    tree: MypyFile | None = None
    # We keep both a list and set of dependencies. A set because it makes it efficient to
    # prevent duplicates and the list because I am afraid of changing the order of
    # iteration over dependencies.
    # They should be managed with add_dependency and suppress_dependency.
    dependencies: list[str]  # Modules directly imported by the module
    dependencies_set: set[str]  # The same but as a set for deduplication purposes
    suppressed: list[str]  # Suppressed/missing dependencies
    suppressed_set: set[str]  # Suppressed/missing dependencies
    priorities: dict[str, int]

    # Map each dependency to the line number where it is first imported
    dep_line_map: dict[str, int]

    # Map from dependency id to its last observed interface hash
    dep_hashes: dict[str, bytes]

    # List of errors reported for this file last time.
    error_lines: list[ErrorTuple]

    # Parent package, its parent, etc.
    ancestors: list[str] | None = None

    # List of (path, line number) tuples giving context for import
    import_context: list[tuple[str, int]]

    # If caller_state is set, the line number in the caller where the import occurred
    caller_line = 0

    # Contains a hash of the public interface in incremental mode
    interface_hash: bytes = b""

    # Hash of import structure that this module depends on. It is not 1:1 with
    # transitive dependencies set, but if two hashes are equal, transitive
    # dependencies are guaranteed to be identical. Some expensive checks can be
    # skipped if this value is unchanged for a module.
    trans_dep_hash: bytes = b""

    # Options, specialized for this file
    options: Options

    # Whether to ignore all errors
    ignore_all = False

    # Errors reported before semantic analysis, to allow fine-grained
    # mode to keep reporting them.
    early_errors: list[ErrorInfo]

    # Type checker used for checking this file.  Use type_checker() for
    # access and to construct this on demand.
    _type_checker: TypeChecker | None = None

    fine_grained_deps_loaded = False

    # Cumulative time spent on this file, in microseconds (for profiling stats)
    time_spent_us: int = 0

    # Per-line type-checking time (cumulative time spent type-checking expressions
    # on a given source code line).
    per_line_checking_time_ns: dict[int, int]

    # Rough estimate of how much time it would take to process this file. Currently,
    # we use file size as a proxy for complexity.
    size_hint: int

    # Mapping from line number to type ignore codes on this line (for imports only).
    imports_ignored: dict[int, list[str]]

    @staticmethod
    def new_state(
        id: str | None,
        path: str | None,
        source: str | None,
        manager: BuildManager,
        caller_state: State | None = None,
        caller_line: int = 0,
        ancestor_for: State | None = None,
        root_source: bool = False,
        # If `temporary` is True, this State is being created to just
        # quickly parse/load the tree, without an intention to further
        # process it. With this flag, any changes to external state as well
        # as error reporting should be avoided.
        temporary: bool = False,
    ) -> State:
        if not temporary:
            assert id or path or source is not None, "Neither id, path nor source given"
        State.order_counter += 1
        if caller_state:
            import_context = caller_state.import_context.copy()
            import_context.append((caller_state.xpath, caller_line))
        else:
            import_context = []
        id = id or "__main__"
        options = manager.options.clone_for_module(id)
        manager.import_options[id] = options.dep_import_options()

        ignore_all = False
        if not path and source is None:
            assert id is not None
            try:
                path, follow_imports = find_module_and_diagnose(
                    manager,
                    id,
                    options,
                    caller_state,
                    caller_line,
                    ancestor_for,
                    root_source,
                    skip_diagnose=temporary,
                )
            except ModuleNotFound as exc:
                if not temporary:
                    manager.missing_modules[id] = exc.reason
                raise
            if follow_imports == "silent":
                ignore_all = True
        elif path and is_silent_import_module(manager, path) and not root_source:
            ignore_all = True

        meta = None
        meta_ex = None
        interface_hash = b""
        meta_source_hash = None
        if path and source is None and manager.cache_enabled:
            meta_pair = find_cache_meta(id, path, manager)
            # TODO: Get mtime if not cached.
            if meta_pair is not None:
                meta, meta_ex = meta_pair
                interface_hash = meta.interface_hash
                meta_source_hash = meta.hash
        if path and source is None and manager.fscache.isdir(path):
            source = ""

        if manager.stats_enabled:
            t0 = time.time()
        meta = validate_meta(meta, id, path, ignore_all, manager)
        if manager.stats_enabled:
            manager.add_stats(validate_meta_time=time.time() - t0)

        if meta:
            assert meta_ex is not None
            # Make copies, since we may modify these and want to
            # compare them to the originals later.
            dependencies = list(meta.dependencies)
            suppressed = list(meta.suppressed)
            all_deps = dependencies + suppressed
            assert len(all_deps) == len(meta.dep_prios)
            priorities = {id: pri for id, pri in zip(all_deps, meta.dep_prios)}
            assert len(all_deps) == len(meta.dep_lines)
            dep_line_map = {id: line for id, line in zip(all_deps, meta.dep_lines)}
            # Merge CacheMetaEx data into CacheMeta data.
            for dep in meta_ex.dependencies + meta_ex.suppressed:
                priorities[dep] = PRI_INDIRECT
                dep_line_map[dep] = 1
            dependencies += meta_ex.dependencies
            meta.dependencies += meta_ex.dependencies
            suppressed += meta_ex.suppressed
            meta.suppressed += meta_ex.suppressed
            meta.dep_hashes += meta_ex.dep_hashes
            assert len(meta.dep_hashes) == len(meta.dependencies)
            dep_hashes = {k: v for (k, v) in zip(meta.dependencies, meta.dep_hashes)}
            # Only copy `error_lines` if the module is not silently imported.
            error_lines = [] if ignore_all else meta_ex.error_lines
            imports_ignored = meta.imports_ignored
        else:
            dependencies = []
            suppressed = []
            priorities = {}
            dep_line_map = {}
            dep_hashes = {}
            error_lines = []
            imports_ignored = {}

        state = State(
            manager=manager,
            order=State.order_counter,
            id=id,
            path=path,
            source=source,
            options=options,
            ignore_all=ignore_all,
            caller_line=caller_line,
            import_context=import_context,
            meta=meta,
            interface_hash=interface_hash,
            meta_source_hash=meta_source_hash,
            dependencies=dependencies,
            suppressed=suppressed,
            priorities=priorities,
            dep_line_map=dep_line_map,
            dep_hashes=dep_hashes,
            error_lines=error_lines,
            imports_ignored=imports_ignored,
        )

        if meta:
            if temporary:
                state.load_tree(temporary=True)
            if not manager.use_fine_grained_cache():
                # Special case: if there were a previously missing package imported here,
                # and it is not present, then we need to re-calculate dependencies.
                # This is to support patterns like this:
                #     from missing_package import missing_module  # type: ignore
                # At first mypy doesn't know that `missing_module` is a module
                # (it may be a variable, a class, or a function), so it is not added to
                # suppressed dependencies. Therefore, when the package with module is added,
                # we need to re-calculate dependencies.
                # NOTE: see comment below for why we skip this in fine-grained mode.
                if exist_added_packages(suppressed, manager):
                    state.needs_parse = True  # This is safe because the cache is anyway stale.
                # This is an inverse to the situation above. If we had an import like this:
                #     from pkg import mod
                # and then mod was deleted, we need to force recompute dependencies, to
                # decide whether we should still depend on a missing pkg.mod. Otherwise,
                # the above import is indistinguishable from something like this:
                #     import pkg
                #     import pkg.mod
                if exist_removed_submodules(dependencies, manager):
                    state.needs_parse = True  # Same as above, the current state is stale anyway.
            state.size_hint = meta.size + MIN_SIZE_HINT
        else:
            # When doing a fine-grained cache load, pretend we only
            # know about modules that have cache information and defer
            # handling new modules until the fine-grained update.
            if manager.use_fine_grained_cache():
                manager.log(f"Deferring module to fine-grained update {path} ({id})")
                raise ModuleNotFound

            if temporary:
                # Eagerly parse temporary states, they are needed rarely.
                state.parse_file(temporary=True)
                state.compute_dependencies()
                if state.manager.workers and state.tree:
                    # We don't need imports in coordinator process anymore, we parse only to
                    # compute dependencies.
                    state.tree.imports = []
                    del state.manager.ast_cache[state.id]
            else:
                state.needs_parse = True

        return state

    def __init__(
        self,
        manager: BuildManager,
        order: int,
        id: str,
        path: str | None,
        source: str | None,
        options: Options,
        ignore_all: bool,
        caller_line: int,
        import_context: list[tuple[str, int]],
        meta: CacheMeta | None,
        interface_hash: bytes,
        meta_source_hash: str | None,
        dependencies: list[str],
        suppressed: list[str],
        priorities: dict[str, int],
        dep_line_map: dict[str, int],
        dep_hashes: dict[str, bytes],
        error_lines: list[ErrorTuple],
        imports_ignored: dict[int, list[str]],
        size_hint: int = 0,
    ) -> None:
        self.manager = manager
        self.order = order
        self.id = id
        self.path = path
        if path:
            # Avoid calling os.abspath, since it makes a getcwd() syscall, which is slow
            if os.path.isabs(path):
                self.abspath = path
            else:
                self.abspath = os.path.normpath(os_path_join(manager.cwd, path))
        self.xpath = path or "<string>"
        self.source = source
        self.options = options
        self.ignore_all = ignore_all
        self.caller_line = caller_line
        self.import_context = import_context
        self.meta = meta
        self.interface_hash = interface_hash
        self.meta_source_hash = meta_source_hash
        self.dependencies = dependencies
        self.suppressed = suppressed
        self.dependencies_set = set(dependencies)
        self.suppressed_set = set(suppressed)
        self.priorities = priorities
        self.dep_line_map = dep_line_map
        self.dep_hashes = dep_hashes
        self.error_lines = error_lines
        self.per_line_checking_time_ns = collections.defaultdict(int)
        self.early_errors = []
        self._type_checker = None
        self.add_ancestors()
        self.imports_ignored = imports_ignored
        self.size_hint = size_hint
        # Pre-computed opaque value of suppressed_deps_opts() used
        # to minimize amount of data sent to parallel workers.
        self.known_suppressed_deps_opts: bytes | None = None
        # An internal flag used by build manager to schedule states for parsing.
        self.needs_parse = False

    def write(self, buf: WriteBuffer) -> None:
        """Serialize State for sending to build worker.

        Note that unlike write() methods for most other classes, this one is
        not idempotent. We erase some bulky values that should either be not needed
        for processing by the worker, or can be re-created from other data relatively
        quickly. These are:
          * self.meta: workers will call self.reload_meta() anyway.
          * self.options: can be restored with Options.clone_for_module().
          * self.error_lines: fresh errors are handled by the coordinator.
        """
        write_int(buf, self.order)
        write_str(buf, self.id)
        write_str_opt(buf, self.path)
        write_str_opt(buf, self.source)  # mostly for mypy -c '<some code>'
        write_bool(buf, self.ignore_all)
        write_int(buf, self.caller_line)
        write_tag(buf, LIST_GEN)
        write_int_bare(buf, len(self.import_context))
        for path, line in self.import_context:
            write_str(buf, path)
            write_int(buf, line)
        write_bytes(buf, self.interface_hash)
        write_str_opt(buf, self.meta_source_hash)
        write_str_list(buf, self.dependencies)
        write_str_list(buf, self.suppressed)
        # TODO: we can possibly serialize these dictionaries in a more compact way.
        # Most keys in the dictionaries should be the same, so we can write them once.
        write_tag(buf, DICT_STR_GEN)
        write_int_bare(buf, len(self.priorities))
        for mod_id, prio in self.priorities.items():
            write_str_bare(buf, mod_id)
            write_int(buf, prio)
        write_tag(buf, DICT_STR_GEN)
        write_int_bare(buf, len(self.dep_line_map))
        for mod_id, line in self.dep_line_map.items():
            write_str_bare(buf, mod_id)
            write_int(buf, line)
        write_tag(buf, DICT_STR_GEN)
        write_int_bare(buf, len(self.dep_hashes))
        for mod_id, dep_hash in self.dep_hashes.items():
            write_str_bare(buf, mod_id)
            write_bytes(buf, dep_hash)
        write_int(buf, self.size_hint)

    @classmethod
    def read(cls, buf: ReadBuffer, manager: BuildManager) -> State:
        order = read_int(buf)
        id = read_str(buf)
        path = read_str_opt(buf)
        source = read_str_opt(buf)
        ignore_all = read_bool(buf)
        caller_line = read_int(buf)
        assert read_tag(buf) == LIST_GEN
        import_context = [(read_str(buf), read_int(buf)) for _ in range(read_int_bare(buf))]
        interface_hash = read_bytes(buf)
        meta_source_hash = read_str_opt(buf)
        dependencies = read_str_list(buf)
        suppressed = read_str_list(buf)
        assert read_tag(buf) == DICT_STR_GEN
        priorities = {read_str_bare(buf): read_int(buf) for _ in range(read_int_bare(buf))}
        assert read_tag(buf) == DICT_STR_GEN
        dep_line_map = {read_str_bare(buf): read_int(buf) for _ in range(read_int_bare(buf))}
        assert read_tag(buf) == DICT_STR_GEN
        dep_hashes = {read_str_bare(buf): read_bytes(buf) for _ in range(read_int_bare(buf))}
        return cls(
            manager=manager,
            order=order,
            id=id,
            path=path,
            source=source,
            # The caller must call clone_for_module().
            options=manager.options,
            ignore_all=ignore_all,
            caller_line=caller_line,
            import_context=import_context,
            meta=None,
            interface_hash=interface_hash,
            meta_source_hash=meta_source_hash,
            dependencies=dependencies,
            suppressed=suppressed,
            priorities=priorities,
            dep_line_map=dep_line_map,
            dep_hashes=dep_hashes,
            error_lines=[],
            imports_ignored={},
            size_hint=read_int(buf),
        )

    def reload_meta(self) -> None:
        """Force reload of cache meta.

        This is used by parallel checking workers to update shared information
        that may have changed after initial graph loading. Currently, this is only
        the interface hash.
        """
        assert self.manager.parallel_worker
        assert self.path is not None
        new_meta_pair = find_cache_meta(self.id, self.path, self.manager, skip_validation=True)
        assert new_meta_pair is not None
        new_meta, _ = new_meta_pair
        # Copy relevant information from new meta (which may be incomplete).
        self.interface_hash = new_meta.interface_hash

    def add_ancestors(self) -> None:
        if self.path is not None:
            _, name = os.path.split(self.path)
            base, _ = os.path.splitext(name)
            if "." in base:
                # This is just a weird filename, don't add anything
                self.ancestors = []
                return
        # All parent packages are new ancestors.
        ancestors = []
        parent = self.id
        while "." in parent:
            parent, _ = parent.rsplit(".", 1)
            ancestors.append(parent)
        self.ancestors = ancestors

    def is_fresh(self) -> bool:
        """Return whether the cache data for this file is fresh."""
        # NOTE: self.dependencies may differ from
        # self.meta.dependencies when a dependency is dropped due to
        # suppression by silent mode.  However, when a suppressed
        # dependency is added back we find out later in the process.
        # Additionally, we need to verify that import following options are
        # same for suppressed dependencies, even if the first check is OK.
        return (
            self.meta is not None
            and self.dependencies == self.meta.dependencies
            and (
                self.options.fine_grained_incremental
                or self.meta.suppressed_deps_opts == self.suppressed_deps_opts()
            )
        )

    def mark_as_rechecked(self) -> None:
        """Marks this module as having been fully re-analyzed by the type-checker."""
        self.manager.rechecked_modules.add(self.id)

    def mark_interface_stale(self) -> None:
        """Marks this module as having a stale public interface, and discards the cache data."""
        self.manager.stale_modules.add(self.id)

    def check_blockers(self) -> None:
        """Raise CompileError if a blocking error is detected."""
        if self.manager.errors.is_blockers():
            self.manager.log("Bailing due to blocking errors")
            self.manager.errors.raise_error()

    @contextlib.contextmanager
    def wrap_context(self, check_blockers: bool = True) -> Iterator[None]:
        """Temporarily change the error import context to match this state.

        Also report an internal error if an unexpected exception was raised
        and raise an exception on a blocking error, unless
        check_blockers is False. Skipping blocking error reporting is used
        in the semantic analyzer so that we can report all blocking errors
        for a file (across multiple targets) to maintain backward
        compatibility.
        """
        save_import_context = self.manager.errors.import_context()
        self.manager.errors.set_import_context(self.import_context)
        try:
            yield
        except CompileError:
            raise
        except Exception as err:
            report_internal_error(
                err,
                self.path,
                0,
                self.manager.errors,
                self.options,
                self.manager.stdout,
                self.manager.stderr,
            )
        self.manager.errors.set_import_context(save_import_context)
        # TODO: Move this away once we've removed the old semantic analyzer?
        if check_blockers:
            self.check_blockers()

    def load_fine_grained_deps(self) -> dict[str, set[str]]:
        return self.manager.load_fine_grained_deps(self.id)

    def load_tree(self, temporary: bool = False) -> None:
        if self.manager.parallel_worker:
            assert self.path is not None
            _, data_file, _ = get_cache_names(self.id, self.path, self.manager.options)
        else:
            assert (
                self.meta is not None
            ), "Internal error: this method must be called only for cached modules"
            data_file = self.meta.data_file

        data: bytes | dict[str, Any] | None
        if self.options.fixed_format_cache:
            data = _load_ff_file(data_file, self.manager, "Could not load tree: ", None)
        else:
            data = _load_json_file(data_file, self.manager, "Load tree ", "Could not load tree: ")
        if data is None:
            return

        t0 = time.time()
        # TODO: Assert data file wasn't changed.
        if isinstance(data, bytes):
            data_io = ReadBuffer(data)
            self.tree = MypyFile.read(data_io)
        else:
            self.tree = MypyFile.deserialize(data)
        t1 = time.time()
        self.manager.add_stats(deserialize_time=t1 - t0)
        if not temporary:
            self.manager.modules[self.id] = self.tree
            self.manager.add_stats(fresh_trees=1)

    def fix_cross_refs(self) -> None:
        assert self.tree is not None, "Internal error: method must be called on parsed file only"
        # Do initial lightweight pass fixing TypeInfos and module cross-references.
        assert modules_state.node_fixer is not None
        modules_state.node_fixer.visit_symbol_table(self.tree.names)
        type_fixer = modules_state.node_fixer.type_fixer
        # Eagerly fix shared instances, before they are used by named_type() calls.
        if instance_cache.str_type is not None:
            instance_cache.str_type.accept(type_fixer)
        if instance_cache.function_type is not None:
            instance_cache.function_type.accept(type_fixer)
        if instance_cache.int_type is not None:
            instance_cache.int_type.accept(type_fixer)
        if instance_cache.bool_type is not None:
            instance_cache.bool_type.accept(type_fixer)
        if instance_cache.object_type is not None:
            instance_cache.object_type.accept(type_fixer)

    # Methods for processing modules from source code.

    def get_source(self) -> str:
        """Get module source and parse inline mypy configurations."""
        manager = self.manager
        t0 = time_ref()

        with self.wrap_context():
            source = self.source
            if self.path and source is None:
                try:
                    path = manager.maybe_swap_for_shadow_path(self.path)
                    source = decode_python_encoding(manager.fscache.read(path))
                    self.source_hash = manager.fscache.hash_digest(path)
                except OSError as ioerr:
                    # ioerr.strerror differs for os.stat failures between Windows and
                    # other systems, but os.strerror(ioerr.errno) does not, so we use that.
                    # (We want the error messages to be platform-independent so that the
                    # tests have predictable output.)
                    assert ioerr.errno is not None
                    raise CompileError(
                        [
                            "mypy: error: Cannot read file '{}': {}".format(
                                self.path.replace(os.getcwd() + os.sep, ""),
                                os.strerror(ioerr.errno),
                            )
                        ],
                        module_with_blocker=self.id,
                    ) from ioerr
                except (UnicodeDecodeError, DecodeError) as decodeerr:
                    if self.path.endswith(".pyd"):
                        err = f"{self.path}: error: Stubgen does not support .pyd files"
                    else:
                        err = f"{self.path}: error: Cannot decode file: {str(decodeerr)}"
                    raise CompileError([err], module_with_blocker=self.id) from decodeerr
            elif self.path and self.manager.fscache.isdir(self.path):
                source = ""
                self.source_hash = ""
            else:
                assert source is not None
                self.source_hash = compute_hash(source)

            self.parse_inline_configuration(source)

            self.size_hint = len(source) + MIN_SIZE_HINT
        self.time_spent_us += time_spent_us(t0)
        return source

    def parse_file_inner(self, source: str, raw_data: FileRawData | None = None) -> None:
        t0 = time_ref()
        self.tree = self.manager.parse_file(
            self.id, self.xpath, source, options=self.options, raw_data=raw_data
        )
        self.time_spent_us += time_spent_us(t0)

    def parse_file(self, *, temporary: bool = False, raw_data: FileRawData | None = None) -> None:
        """Parse file and run first pass of semantic analysis.

        Everything done here is local to the file. Don't depend on imported
        modules in any way. Logic here should be kept in sync with BuildManager.parse_all().
        """
        self.needs_parse = False
        tree = self.tree
        if tree is not None:
            # The file was already parsed.
            return

        if raw_data is None:
            source = self.get_source()
        else:
            source = ""
        manager = self.manager
        # Can we reuse a previously parsed AST? This avoids redundant work in daemon.
        if self.id not in manager.ast_cache:
            self.manager.log(f"Parsing {self.xpath} ({self.id})")
            ignore_errors = self.ignore_all or self.options.ignore_errors
            if ignore_errors:
                self.manager.errors.ignored_files.add(self.xpath)
            with self.wrap_context():
                manager.errors.set_file(self.xpath, self.id, options=self.options)
                if raw_data is not None:
                    # Apply inline mypy config before deserialization, since
                    # some options (e.g. implicit_optional) affect how the
                    # AST is built during deserialization.
                    self.source_hash = raw_data.source_hash
                    self.apply_inline_configuration(raw_data.mypy_comments)
                self.parse_file_inner(source, raw_data)
                assert self.tree is not None
                # New parser returns serialized trees that need to be de-serialized.
                if self.tree.raw_data is not None:
                    assert raw_data is None
                    self.tree = load_from_raw(
                        self.xpath,
                        self.id,
                        self.tree.raw_data,
                        manager.errors,
                        self.options,
                        imports_only=bool(self.manager.workers),
                    )
                if manager.errors.is_blockers():
                    manager.log("Bailing due to parse errors")
                    manager.errors.raise_error()
            # Make a copy of any errors produced during parse time so that
            # fine-grained mode can repeat them when the module is
            # reprocessed.
            self.early_errors = list(manager.errors.error_info_map.get(self.xpath, []))
            self.semantic_analysis_pass1()
        else:
            # Reuse a cached AST
            manager.log(f"Using cached AST for {self.xpath} ({self.id})")
            self.tree, self.early_errors, source_hash = manager.ast_cache[self.id]
            self.source_hash = source_hash

        assert self.tree is not None
        if not temporary:
            manager.modules[self.id] = self.tree
            self.check_blockers()

        manager.ast_cache[self.id] = (self.tree, self.early_errors, self.source_hash)
        assert self.tree is not None
        if self.tree.raw_data is not None:
            # Size of serialized tree is a better proxy for file complexity than
            # file size, so we use that when possible. Note that we rely on lucky
            # coincidence that serialized tree size has same order of magnitude as
            # file size, so we don't need any normalization factor in situations
            # where parsed and cached files are mixed.
            self.size_hint = len(self.tree.raw_data.defs) + MIN_SIZE_HINT
        self.setup_errors()

    def setup_errors(self) -> None:
        assert self.tree is not None
        self.manager.errors.set_file_ignored_lines(
            self.xpath, self.tree.ignored_lines, self.ignore_all or self.options.ignore_errors
        )
        self.manager.errors.set_skipped_lines(self.xpath, self.tree.skipped_lines)

    def parse_inline_configuration(self, source: str) -> None:
        """Check for inline mypy: options directive and parse them."""
        flags = get_mypy_comments(source)
        self.apply_inline_configuration(flags)

    def apply_inline_configuration(self, flags: list[tuple[int, str]] | None) -> None:
        """Apply inline mypy configuration comments and check for invalid options."""
        if flags:
            changes, config_errors = parse_mypy_comments(flags, self.options)
            self.options = self.options.apply_changes(changes)
            self.manager.errors.set_file(self.xpath, self.id, self.options)
            for lineno, error in config_errors:
                self.manager.error(lineno, error)
        self.check_for_invalid_options()

    def check_for_invalid_options(self) -> None:
        if self.options.mypyc and not self.options.strict_bytes:
            self.manager.errors.set_file(self.xpath, self.id, options=self.options)
            self.manager.error(
                None, "Option --strict-bytes cannot be disabled when using mypyc", blocker=True
            )

    def semantic_analysis_pass1(self) -> None:
        """Perform pass 1 of semantic analysis, which happens immediately after parsing.

        This pass can't assume that any other modules have been processed yet.
        """
        options = self.options
        assert self.tree is not None

        t0 = time_ref()

        # Do the first pass of semantic analysis: analyze the reachability
        # of blocks and import statements. We must do this before
        # processing imports, since this may mark some import statements as
        # unreachable.
        #
        # TODO: This should not be considered as a semantic analysis
        #     pass -- it's an independent pass.
        if not options.native_parser or not self.manager.fscache.exists(
            self.xpath, real_only=True
        ):
            analyzer = SemanticAnalyzerPreAnalysis()
            with self.wrap_context():
                analyzer.visit_file(self.tree, self.xpath, self.id, options)
        # TODO: Do this while constructing the AST?
        self.tree.names = SymbolTable()
        if not self.tree.is_stub:
            if not self.options.allow_redefinition:
                # Perform some low-key variable renaming when assignments can't
                # widen inferred types
                self.tree.accept(LimitedVariableRenameVisitor())
            if options.allow_redefinition_old:
                # Perform more renaming across the AST to allow variable redefinitions
                self.tree.accept(VariableRenameVisitor())
        self.time_spent_us += time_spent_us(t0)

    def add_dependency(self, dep: str) -> None:
        if dep not in self.dependencies_set:
            self.dependencies.append(dep)
            self.dependencies_set.add(dep)
        if dep in self.suppressed_set:
            self.suppressed.remove(dep)
            self.suppressed_set.remove(dep)

    def suppress_dependency(self, dep: str) -> None:
        if dep in self.dependencies_set:
            self.dependencies.remove(dep)
            self.dependencies_set.remove(dep)
        if dep not in self.suppressed_set:
            self.suppressed.append(dep)
            self.suppressed_set.add(dep)

    def compute_dependencies(self) -> None:
        """Compute a module's dependencies after parsing it.

        This is used when we parse a file that we didn't have
        up-to-date cache information for. When we have an up-to-date
        cache, we just use the cached info.
        """
        manager = self.manager
        assert self.tree is not None

        # Compute (direct) dependencies.
        # Add all direct imports (this is why we needed the first pass).
        # Also keep track of each dependency's source line.
        # Missing dependencies will be moved from dependencies to
        # suppressed when they fail to be loaded in load_graph.

        self.dependencies = []
        self.dependencies_set = set()
        self.suppressed = []
        self.suppressed_set = set()
        self.priorities = {}  # id -> priority
        self.dep_line_map = {}  # id -> line
        self.dep_hashes = {}
        # We copy imports as defs to (partially) support some legacy mypy plugins,
        # most notably old NumPy plugin that does some imports patching, see #21323.
        copied_imports = False
        if not self.tree.defs and self.tree.raw_data is not None:
            self.tree.defs = list(self.tree.imports)
            copied_imports = True
        dep_entries = manager.all_imported_modules_in_file(
            self.tree
        ) + self.manager.plugin.get_additional_deps(self.tree)
        if copied_imports:
            self.tree.defs = []
        for pri, id, line in dep_entries:
            self.priorities[id] = min(pri, self.priorities.get(id, PRI_ALL))
            if id == self.id:
                continue
            self.add_dependency(id)
            if id not in self.dep_line_map:
                self.dep_line_map[id] = line
        import_lines = self.dep_line_map.values()
        self.imports_ignored = {
            line: codes for line, codes in self.tree.ignored_lines.items() if line in import_lines
        }
        # Every module implicitly depends on builtins.
        if self.id != "builtins":
            self.add_dependency("builtins")
        if self.tree.uses_template_strings:
            self.add_dependency("string.templatelib")

        self.check_blockers()  # Can fail due to bogus relative imports

    def type_check_first_pass(self, recurse_into_functions: bool = True) -> None:
        if self.options.semantic_analysis_only:
            return
        t0 = time_ref()
        with self.wrap_context():
            self.type_checker().check_first_pass(recurse_into_functions=recurse_into_functions)
        self.time_spent_us += time_spent_us(t0)

    def type_checker(self) -> TypeChecker:
        if not self._type_checker:
            assert self.tree is not None, "Internal error: must be called on parsed file only"
            manager = self.manager
            self._type_checker = TypeChecker(
                manager.errors,
                manager.modules,
                self.options,
                self.tree,
                self.xpath,
                manager.plugin,
                self.per_line_checking_time_ns,
            )
        return self._type_checker

    def type_map(self) -> dict[Expression, Type]:
        # We can extract the master type map directly since at this
        # point no temporary type maps can be active.
        assert len(self.type_checker()._type_maps) == 1
        return self.type_checker()._type_maps[0]

    def type_check_second_pass(
        self,
        todo: Sequence[DeferredNode] | None = None,
        recurse_into_functions: bool = True,
        impl_only: bool = False,
    ) -> bool:
        if self.options.semantic_analysis_only:
            return False
        t0 = time_ref()
        with self.wrap_context():
            result = self.type_checker().check_second_pass(
                todo=todo, recurse_into_functions=recurse_into_functions, impl_only=impl_only
            )
        self.time_spent_us += time_spent_us(t0)
        return result

    def detect_possibly_undefined_vars(self) -> None:
        assert self.tree is not None, "Internal error: method must be called on parsed file only"
        if self.tree.is_stub:
            # We skip stub files because they aren't actually executed.
            return
        manager = self.manager
        manager.errors.set_file(self.xpath, self.tree.fullname, options=self.options)
        if manager.errors.is_error_code_enabled(
            codes.POSSIBLY_UNDEFINED
        ) or manager.errors.is_error_code_enabled(codes.USED_BEFORE_DEF):
            with self.wrap_context():
                self.tree.accept(
                    PossiblyUndefinedVariableVisitor(
                        MessageBuilder(manager.errors, manager.modules),
                        self.type_map(),
                        self.options,
                        self.tree.names,
                    )
                )

    def finish_passes(self) -> None:
        assert self.tree is not None, "Internal error: method must be called on parsed file only"
        manager = self.manager
        if self.options.semantic_analysis_only:
            return
        t0 = time_ref()
        with self.wrap_context():
            # Some tests (and tools) want to look at the set of all types.
            options = manager.options
            if options.export_types:
                manager.all_types.update(self.type_map())

            # We should always patch indirect dependencies, even in full (non-incremental) builds,
            # because the cache still may be written, and it must be correct.
            self.patch_indirect_dependencies(
                # Two possible sources of indirect dependencies:
                # * Symbols not directly imported in this module but accessed via an attribute
                #   or via a re-export (vast majority of these recorded in semantic analysis).
                # * For each expression type we need to record definitions of type components
                #   since "meaning" of the type may be updated when definitions are updated.
                self.tree.module_refs | self.type_checker().module_refs,
                set(self.type_map().values()),
            )

            if self.options.dump_inference_stats:
                dump_type_stats(
                    self.tree,
                    self.xpath,
                    modules=self.manager.modules,
                    inferred=True,
                    typemap=self.type_map(),
                )
            manager.report_file(self.tree, self.type_map(), self.options)

            self.update_fine_grained_deps(self.manager.fg_deps)

            if manager.options.export_ref_info:
                write_undocumented_ref_info(
                    self, manager.metastore, manager.options, self.type_map()
                )

            self.free_state()
            if not manager.options.fine_grained_incremental and not manager.options.preserve_asts:
                free_tree(self.tree)
                self.tree.defs.clear()
        self.time_spent_us += time_spent_us(t0)

    def free_state(self) -> None:
        if self._type_checker:
            self._type_checker.reset()
            self._type_checker = None

    def patch_indirect_dependencies(self, module_refs: set[str], types: set[Type]) -> None:
        assert self.ancestors is not None
        existing_deps = set(self.dependencies + self.suppressed + self.ancestors)
        existing_deps.add(self.id)

        encountered = self.manager.indirection_detector.find_modules(types) | module_refs
        for dep in sorted(encountered - existing_deps):
            if dep not in self.manager.modules:
                continue
            self.add_dependency(dep)
            self.priorities[dep] = PRI_INDIRECT

    def compute_fine_grained_deps(self) -> dict[str, set[str]]:
        assert self.tree is not None
        if self.id in ("builtins", "typing", "types", "sys", "_typeshed"):
            # We don't track changes to core parts of typeshed -- the
            # assumption is that they are only changed as part of mypy
            # updates, which will invalidate everything anyway. These
            # will always be processed in the initial non-fine-grained
            # build. Other modules may be brought in as a result of an
            # fine-grained increment, and we may need these
            # dependencies then to handle cyclic imports.
            return {}
        from mypy.server.deps import get_dependencies  # Lazy import to speed up startup

        return get_dependencies(
            target=self.tree,
            type_map=self.type_map(),
            python_version=self.options.python_version,
            options=self.manager.options,
        )

    def update_fine_grained_deps(self, deps: dict[str, set[str]]) -> None:
        options = self.manager.options
        if options.cache_fine_grained or options.fine_grained_incremental:
            from mypy.server.deps import merge_dependencies  # Lazy import to speed up startup

            merge_dependencies(self.compute_fine_grained_deps(), deps)
            type_state.update_protocol_deps(deps)

    def suppressed_deps_opts(self) -> bytes:
        if not self.suppressed:
            return b""
        if self.known_suppressed_deps_opts:
            return self.known_suppressed_deps_opts
        buf = WriteBuffer()
        import_options = self.manager.import_options
        for dep in sorted(self.suppressed):
            # Using .get() is a bit defensive, but just in case we have a bug elsewhere
            # (e.g. in the daemon), it is better to get a stale cache than a crash.
            reason = self.manager.missing_modules.get(dep, SuppressionReason.NOT_FOUND)
            if self.priorities.get(dep) != PRI_INDIRECT:
                write_str_bare(buf, dep)
                write_bytes_bare(buf, import_options[dep])
                write_int_bare(buf, reason)
        return buf.getvalue()

    def write_cache(self) -> tuple[CacheMeta, str] | None:
        assert self.tree is not None, "Internal error: method must be called on parsed file only"
        # We don't support writing cache files in fine-grained incremental mode.
        if (
            not self.path
            or self.options.cache_dir == os.devnull
            or self.options.fine_grained_incremental
        ):
            if self.options.debug_serialize:
                try:
                    if self.manager.options.fixed_format_cache:
                        data = WriteBuffer()
                        self.tree.write(data)
                    else:
                        self.tree.serialize()
                except Exception:
                    print(f"Error serializing {self.id}", file=self.manager.stdout)
                    raise  # Propagate to display traceback
            return None
        dep_prios = self.dependency_priorities()
        dep_lines = self.dependency_lines()
        assert self.source_hash is not None
        assert len(set(self.dependencies)) == len(
            self.dependencies
        ), f"Duplicates in dependencies list for {self.id} ({self.dependencies})"
        new_interface_hash, meta_tuple = write_cache(
            self.id,
            self.path,
            self.tree,
            # Indirect dependencies are stored separately as part of CacheMetaEx.
            [dep for dep in self.dependencies if self.priorities.get(dep) != PRI_INDIRECT],
            [dep for dep in self.suppressed if self.priorities.get(dep) != PRI_INDIRECT],
            self.suppressed_deps_opts(),
            self.imports_ignored,
            dep_prios,
            dep_lines,
            self.interface_hash,
            self.trans_dep_hash,
            self.source_hash,
            self.ignore_all,
            self.manager,
        )
        if new_interface_hash == self.interface_hash:
            self.manager.log(f"Cached module {self.id} has same interface")
        else:
            self.manager.log(f"Cached module {self.id} has changed interface")
            self.mark_interface_stale()
            self.interface_hash = new_interface_hash
        return meta_tuple

    def verify_dependencies(self, suppressed_only: bool = False) -> None:
        """Report errors for import targets in modules that don't exist.

        If suppressed_only is set, only check suppressed dependencies.
        """
        manager = self.manager
        assert self.ancestors is not None
        # Strip out indirect dependencies. See comment in build.load_graph().
        if suppressed_only:
            all_deps = [dep for dep in self.suppressed if self.priorities.get(dep) != PRI_INDIRECT]
        else:
            dependencies = [
                dep
                for dep in self.dependencies + self.suppressed
                if self.priorities.get(dep) != PRI_INDIRECT
            ]
            all_deps = dependencies + self.ancestors
        for dep in all_deps:
            if dep in manager.modules:
                continue
            options = manager.options.clone_for_module(dep)
            if options.ignore_missing_imports:
                continue
            line = self.dep_line_map.get(dep, 1)
            try:
                if dep in self.ancestors:
                    state: State | None = None
                    ancestor: State | None = self
                else:
                    state, ancestor = self, None
                # Called just for its side effects of producing diagnostics.
                find_module_and_diagnose(
                    manager,
                    dep,
                    options,
                    caller_state=state,
                    caller_line=line,
                    ancestor_for=ancestor,
                )
            except (ModuleNotFound, CompileError):
                # Swallow up any ModuleNotFounds or CompilerErrors while generating
                # a diagnostic. CompileErrors may get generated in
                # fine-grained mode when an __init__.py is deleted, if a module
                # that was in that package has targets reprocessed before
                # it is renamed.
                pass

    def dependency_priorities(self) -> list[int]:
        return [
            prio
            for dep in self.dependencies + self.suppressed
            if (prio := self.priorities.get(dep, PRI_HIGH)) != PRI_INDIRECT
        ]

    def dependency_lines(self) -> list[int]:
        return [
            self.dep_line_map.get(dep, 1)
            for dep in self.dependencies + self.suppressed
            if self.priorities.get(dep) != PRI_INDIRECT
        ]

    def generate_unused_ignore_notes(self) -> None:
        if (
            self.options.warn_unused_ignores
            or codes.UNUSED_IGNORE in self.options.enabled_error_codes
        ) and codes.UNUSED_IGNORE not in self.options.disabled_error_codes:
            # We only need this for the daemon, regular incremental does this unconditionally.
            if self.meta and self.options.fine_grained_incremental:
                self.verify_dependencies(suppressed_only=True)
            is_typeshed = self.tree is not None and self.tree.is_typeshed_file(self.options)
            with self.wrap_context():
                self.manager.errors.generate_unused_ignore_errors(self.xpath, is_typeshed)

    def generate_ignore_without_code_notes(self) -> None:
        if self.manager.errors.is_error_code_enabled(codes.IGNORE_WITHOUT_CODE):
            is_typeshed = self.tree is not None and self.tree.is_typeshed_file(self.options)
            with self.wrap_context():
                self.manager.errors.generate_ignore_without_code_errors(
                    self.xpath, self.options.warn_unused_ignores, is_typeshed
                )


# Module import and diagnostic glue


def find_module_and_diagnose(
    manager: BuildManager,
    id: str,
    options: Options,
    caller_state: State | None = None,
    caller_line: int = 0,
    ancestor_for: State | None = None,
    root_source: bool = False,
    skip_diagnose: bool = False,
) -> tuple[str, str]:
    """Find a module by name, respecting follow_imports and producing diagnostics.

    If the module is not found, then the ModuleNotFound exception is raised.

    Args:
      id: module to find
      options: the options for the module being loaded
      caller_state: the state of the importing module, if applicable
      caller_line: the line number of the import
      ancestor_for: the child module this is an ancestor of, if applicable
      root_source: whether this source was specified on the command line
      skip_diagnose: skip any error diagnosis and reporting (but ModuleNotFound is
          still raised if the module is missing)

    The specified value of follow_imports for a module can be overridden
    if the module is specified on the command line or if it is a stub,
    so we compute and return the "effective" follow_imports of the module.

    Returns a tuple containing (file path, target's effective follow_imports setting)
    """
    result = find_module_with_reason(id, manager)
    if isinstance(result, str):
        # For non-stubs, look at options.follow_imports:
        # - normal (default) -> fully analyze
        # - silent -> analyze but silence errors
        # - skip -> don't analyze, make the type Any
        follow_imports = options.follow_imports
        if (
            root_source  # Honor top-level modules
            or (
                result.endswith(".pyi")  # Stubs are always normal
                and not options.follow_imports_for_stubs  # except when they aren't
            )
            or id in CORE_BUILTIN_MODULES  # core is always normal
        ):
            follow_imports = "normal"
        if skip_diagnose:
            pass
        elif follow_imports == "silent":
            # Still import it, but silence non-blocker errors.
            manager.log(f"Silencing {result} ({id})")
        elif follow_imports == "skip" or follow_imports == "error":
            # In 'error' mode, produce special error messages.
            if id not in manager.missing_modules:
                manager.log(f"Skipping {result} ({id})")
            if follow_imports == "error":
                if ancestor_for:
                    skipping_ancestor(manager, id, result, ancestor_for)
                else:
                    skipping_module(manager, caller_line, caller_state, id, result)
            reason = SuppressionReason.SKIPPED
            if options.ignore_missing_imports:
                # Performance optimization: when we are ignoring imports, there is no
                # difference for the caller between skipped import and actually missing one.
                reason = SuppressionReason.NOT_FOUND
            raise ModuleNotFound(reason=reason)
        if is_silent_import_module(manager, result) and not root_source:
            follow_imports = "silent"
        return result, follow_imports
    else:
        # Could not find a module.  Typically, the reason is a
        # misspelled module name, missing stub, module not in
        # search path or the module has not been installed.

        ignore_missing_imports = options.ignore_missing_imports

        if skip_diagnose:
            raise ModuleNotFound
        if caller_state:
            if not (ignore_missing_imports or in_partial_package(id, manager)):
                module_not_found(manager, caller_line, caller_state, id, result)
            raise ModuleNotFound
        elif root_source:
            # If we can't find a root source it's always fatal.
            # TODO: This might hide non-fatal errors from
            # root sources processed earlier.
            raise CompileError([f'mypy: error: Cannot find module "{id}"'])
        else:
            raise ModuleNotFound


def exist_added_packages(suppressed: list[str], manager: BuildManager) -> bool:
    """Find if there are any newly added packages that were previously suppressed.

    Exclude everything not in build for follow-imports=skip.
    """
    for dep in suppressed:
        if dep in manager.source_set.source_modules:
            # We don't need to add any special logic for this. If a module
            # is added to build, importers will be invalidated by normal mechanism.
            continue
        path = find_module_simple(dep, manager)
        if not path:
            continue
        options = manager.options.clone_for_module(dep)
        # Technically this is not 100% correct, since we can have:
        #     from pkg import mod
        # with
        #     [mypy-pkg]
        #     follow-import = silent
        #     [mypy-pkg.mod]
        #     follow-imports = normal
        # But such cases are extremely rare, and this allows us to avoid
        # massive performance impact in much more common situations.
        if options.follow_imports in ("skip", "error") and (
            not path.endswith(".pyi") or options.follow_imports_for_stubs
        ):
            continue
        if os.path.basename(path) in ("__init__.py", "__init__.pyi"):
            return True
    return False


def exist_removed_submodules(dependencies: list[str], manager: BuildManager) -> bool:
    """Find if there are any submodules of packages that are now missing.

    This is conceptually an inverse of exist_added_packages().
    """
    dependencies_set = set(dependencies)
    for dep in dependencies:
        if "." not in dep:
            continue
        if dep in manager.source_set.source_modules:
            # We still know it is definitely a module.
            continue
        direct_ancestor, _ = dep.rsplit(".", maxsplit=1)
        if direct_ancestor not in dependencies_set:
            continue
        if find_module_simple(dep, manager) is None:
            return True
    return False


def find_module_simple(id: str, manager: BuildManager) -> str | None:
    """Find a filesystem path for module `id` or `None` if not found."""
    if manager.stats_enabled:
        t0 = time.time()
    x = manager.find_module_cache.find_module(id, fast_path=True)
    if manager.stats_enabled:
        manager.add_stats(find_module_time=time.time() - t0, find_module_calls=1)
    if isinstance(x, ModuleNotFoundReason):
        return None
    return x


def find_module_with_reason(id: str, manager: BuildManager) -> ModuleSearchResult:
    """Find a filesystem path for module `id` or the reason it can't be found."""
    if manager.stats_enabled:
        t0 = time.time()
    x = manager.find_module_cache.find_module(id, fast_path=False)
    if manager.stats_enabled:
        manager.add_stats(find_module_time=time.time() - t0, find_module_calls=1)
    return x


def in_partial_package(id: str, manager: BuildManager) -> bool:
    """Check if a missing module can potentially be a part of a package.

    This checks if there is any existing parent __init__.pyi stub that
    defines a module-level __getattr__ (a.k.a. partial stub package).
    """
    while "." in id:
        ancestor, _ = id.rsplit(".", 1)
        if ancestor in manager.known_partial_packages:
            return manager.known_partial_packages[ancestor]
        if ancestor in manager.modules:
            ancestor_mod: MypyFile | None = manager.modules[ancestor]
        else:
            # Ancestor is not in build, try quickly if we can find it.
            try:
                ancestor_st = State.new_state(
                    id=ancestor, path=None, source=None, manager=manager, temporary=True
                )
            except (ModuleNotFound, CompileError):
                ancestor_mod = None
            else:
                ancestor_mod = ancestor_st.tree
                # We will not need this anymore.
                ancestor_st.tree = None
        if ancestor_mod is not None:
            # Bail out soon, complete subpackage found
            manager.known_partial_packages[ancestor] = ancestor_mod.is_partial_stub_package
            return ancestor_mod.is_partial_stub_package
        id = ancestor
    return False


def module_not_found(
    manager: BuildManager,
    line: int,
    caller_state: State,
    target: str,
    reason: ModuleNotFoundReason,
) -> None:
    errors = manager.errors
    save_import_context = errors.import_context()
    errors.set_import_context(caller_state.import_context)
    errors.set_file(caller_state.xpath, caller_state.id, caller_state.options)
    errors.set_file_ignored_lines(
        caller_state.xpath,
        caller_state.tree.ignored_lines if caller_state.tree else caller_state.imports_ignored,
        caller_state.ignore_all or caller_state.options.ignore_errors,
    )
    if target == "builtins":
        manager.error(
            line, 'Cannot find "builtins" module. Typeshed appears broken!', blocker=True
        )
        errors.raise_error()
    else:
        daemon = manager.options.fine_grained_incremental
        msg, notes = reason.error_message_templates(daemon)
        if reason == ModuleNotFoundReason.NOT_FOUND:
            code = codes.IMPORT_NOT_FOUND
        elif (
            reason == ModuleNotFoundReason.FOUND_WITHOUT_TYPE_HINTS
            or reason == ModuleNotFoundReason.APPROVED_STUBS_NOT_INSTALLED
        ):
            code = codes.IMPORT_UNTYPED
        else:
            code = codes.IMPORT
        manager.error(line, msg.format(module=target), code=code)

        if (
            reason == ModuleNotFoundReason.NOT_FOUND
            and not errors.prefer_simple_messages()
            and errors.is_error_code_enabled(code)
            and line not in errors.ignored_lines.get(caller_state.xpath, {})
        ):
            top_level_target = target.split(".")[0]
            if not top_level_target.startswith("_"):
                known_modules = get_known_modules(
                    manager.find_module_cache.stdlib_py_versions, manager.options.python_version
                )
                matches = best_matches(top_level_target, known_modules, n=3)
                matches = [m for m in matches if m.lower() != top_level_target.lower()]
                if matches:
                    errors.report(
                        line,
                        0,
                        f'Did you mean {pretty_seq(matches, "or")}?',
                        severity="note",
                        code=code,
                    )

        dist = stub_distribution_name(target)
        for note in notes:
            if "{stub_dist}" in note:
                assert dist is not None
                note = note.format(stub_dist=dist)
            manager.note(line, note, only_once=True, code=code)
        if reason is ModuleNotFoundReason.APPROVED_STUBS_NOT_INSTALLED:
            assert dist is not None
            manager.missing_stub_packages.add(dist)
    errors.set_import_context(save_import_context)


def skipping_module(
    manager: BuildManager, line: int, caller_state: State | None, id: str, path: str
) -> None:
    """Produce an error for an import ignored due to --follow_imports=error"""
    assert caller_state, (id, path)
    save_import_context = manager.errors.import_context()
    manager.errors.set_import_context(caller_state.import_context)
    manager.errors.set_file(caller_state.xpath, caller_state.id, manager.options)
    manager.error(line, f'Import of "{id}" ignored')
    manager.note(
        line, "(Using --follow-imports=error, module not passed on command line)", only_once=True
    )
    manager.errors.set_import_context(save_import_context)


def skipping_ancestor(manager: BuildManager, id: str, path: str, ancestor_for: State) -> None:
    """Produce an error for an ancestor ignored due to --follow_imports=error"""
    # TODO: Read the path (the __init__.py file) and return
    # immediately if it's empty or only contains comments.
    # But beware, some package may be the ancestor of many modules,
    # so we'd need to cache the decision.
    save_import_context = manager.errors.import_context()
    manager.errors.set_import_context([])
    manager.errors.set_file(ancestor_for.xpath, ancestor_for.id, manager.options)
    manager.error(None, f'Ancestor package "{id}" ignored', only_once=True)
    manager.note(
        None, "(Using --follow-imports=error, submodule passed on command line)", only_once=True
    )
    manager.errors.set_import_context(save_import_context)


def log_configuration(manager: BuildManager, sources: list[BuildSource]) -> None:
    """Output useful configuration information to LOG and TRACE"""

    if not manager.logging_enabled:
        return

    config_file = manager.options.config_file
    if config_file:
        config_file = os.path.abspath(config_file)

    manager.log()
    configuration_vars = [
        ("Mypy Version", __version__),
        ("Config File", (config_file or "Default")),
        ("Configured Executable", manager.options.python_executable or "None"),
        ("Current Executable", sys.executable),
        ("Cache Dir", manager.options.cache_dir),
        ("Compiled", str(not __file__.endswith(".py"))),
        ("Exclude", manager.options.exclude),
    ]

    for conf_name, conf_value in configuration_vars:
        manager.log(f"{conf_name + ':':24}{conf_value}")

    for source in sources:
        manager.log(f"{'Found source:':24}{source}")

    # Complete list of searched paths can get very long, put them under TRACE
    for path_type, paths in manager.search_paths.asdict().items():
        if not paths:
            manager.trace(f"No {path_type}")
            continue

        manager.trace(f"{path_type}:")

        for pth in paths:
            manager.trace(f"    {pth}")


# The driver


def dispatch(
    sources: list[BuildSource],
    manager: BuildManager,
    stdout: TextIO,
    connect_threads: list[Thread],
) -> Graph:
    log_configuration(manager, sources)

    t0 = time.time()

    # We disable GC while loading the graph as a performance optimization for
    # cold-cache runs. The parsed ASTs are trees, and therefore should not have any
    # reference cycles. This is an important optimization, since we create a lot of
    # new objects while parsing files.
    global initial_gc_freeze_done
    if (
        not manager.options.test_env
        and platform.python_implementation() == "CPython"
        and not initial_gc_freeze_done
    ):
        gc.disable()
    graph = load_graph(sources, manager)

    # This is a kind of unfortunate hack to work around some of fine-grained's
    # fragility: if we have loaded less than 50% of the specified files from
    # cache in fine-grained cache mode, load the graph again honestly.
    # In this case, we just turn the cache off entirely, so we don't need
    # to worry about some files being loaded and some from cache and so
    # that fine-grained mode never *writes* to the cache.
    if manager.use_fine_grained_cache() and len(graph) < 0.50 * len(sources):
        manager.log("Redoing load_graph without cache because too much was missing")
        manager.cache_enabled = False
        graph = load_graph(sources, manager)

    if (
        not manager.options.test_env
        and platform.python_implementation() == "CPython"
        and not initial_gc_freeze_done
    ):
        gc.freeze()
        gc.unfreeze()
        gc.enable()
        initial_gc_freeze_done = True

    for id in graph:
        manager.import_map[id] = graph[id].dependencies_set

    t1 = time.time()
    manager.add_stats(
        graph_size=len(graph),
        stubs_found=sum(g.path is not None and g.path.endswith(".pyi") for g in graph.values()),
        graph_load_time=(t1 - t0),
        fm_cache_size=len(manager.find_module_cache.results),
    )
    if not graph:
        print("Nothing to do?!", file=stdout)
        return graph
    manager.log(f"Loaded graph with {len(graph)} nodes ({t1 - t0:.3f} sec)")
    if manager.options.dump_graph:
        dump_graph(graph, stdout)
        return graph

    # Fine-grained dependencies that didn't have an associated module in the build
    # are serialized separately, so we read them after we load the graph.
    # We need to read them both for running in daemon mode and if we are generating
    # a fine-grained cache (so that we can properly update them incrementally).
    # The `read_deps_cache` will also validate
    # the deps cache against the loaded individual cache files.
    if manager.options.cache_fine_grained or manager.use_fine_grained_cache():
        t2 = time.time()
        fg_deps_meta = read_deps_cache(manager, graph)
        manager.add_stats(load_fg_deps_time=time.time() - t2)
        if fg_deps_meta is not None:
            manager.fg_deps_meta = fg_deps_meta
        elif manager.stats.get("fresh_metas", 0) > 0:
            # Clear the stats, so we don't infinite loop because of positive fresh_metas
            manager.stats.clear()
            # There were some cache files read, but no fine-grained dependencies loaded.
            manager.log("Error reading fine-grained dependencies cache -- aborting cache load")
            manager.cache_enabled = False
            manager.log("Falling back to full run -- reloading graph...")
            return dispatch(sources, manager, stdout, connect_threads)

    # If we are loading a fine-grained incremental mode cache, we
    # don't want to do a real incremental reprocess of the
    # graph---we'll handle it all later.
    if not manager.use_fine_grained_cache():
        # Wait for workers since they may be needed at this point.
        for thread in connect_threads:
            thread.join()
        not_connected = [str(idx) for idx, wc in enumerate(manager.workers) if not wc.connected]
        if not_connected:
            raise OSError(f"Cannot connect to build worker(s): {', '.join(not_connected)}")
        process_graph(graph, manager)
        # Update plugins snapshot.
        write_plugins_snapshot(manager)
        manager.old_plugins_snapshot = manager.plugins_snapshot
        if manager.options.cache_fine_grained or manager.options.fine_grained_incremental:
            # If we are running a daemon or are going to write cache for further fine-grained use,
            # then we need to collect fine-grained protocol dependencies.
            # Since these are a global property of the program, they are calculated after we
            # processed the whole graph.
            type_state.add_all_protocol_deps(manager.fg_deps)
            if not manager.options.fine_grained_incremental:
                rdeps = generate_deps_for_cache(manager, graph)
                write_deps_cache(rdeps, manager, graph)

    if manager.options.dump_deps:
        # This speeds up startup a little when not using the daemon mode.
        from mypy.server.deps import dump_all_dependencies

        dump_all_dependencies(
            manager.modules, manager.all_types, manager.options.python_version, manager.options
        )

    return graph


class NodeInfo:
    """Some info about a node in the graph of SCCs."""

    def __init__(self, index: int, scc: list[str]) -> None:
        self.node_id = "n%d" % index
        self.scc = scc
        self.sizes: dict[str, int] = {}  # mod -> size in bytes
        self.deps: dict[str, int] = {}  # node_id -> pri

    def dumps(self) -> str:
        """Convert to JSON string."""
        total_size = sum(self.sizes.values())
        return "[{}, {}, {},\n     {},\n     {}]".format(
            json.dumps(self.node_id),
            json.dumps(total_size),
            json.dumps(self.scc),
            json.dumps(self.sizes),
            json.dumps(self.deps),
        )


def dump_timing_stats(path: str, graph: Graph) -> None:
    """Dump timing stats for each file in the given graph."""
    with open(path, "w") as f:
        for id in sorted(graph):
            f.write(f"{id} {graph[id].time_spent_us}\n")


def dump_line_checking_stats(path: str, graph: Graph) -> None:
    """Dump per-line expression type checking stats."""
    with open(path, "w") as f:
        for id in sorted(graph):
            if not graph[id].per_line_checking_time_ns:
                continue
            f.write(f"{id}:\n")
            for line in sorted(graph[id].per_line_checking_time_ns):
                line_time = graph[id].per_line_checking_time_ns[line]
                f.write(f"{line:>5} {line_time/1000:8.1f}\n")


def dump_graph(graph: Graph, stdout: TextIO | None = None) -> None:
    """Dump the graph as a JSON string to stdout.

    This copies some of the work by process_graph()
    (sorted_components() and order_ascc()).
    """
    stdout = stdout or sys.stdout
    nodes = []
    sccs = sorted_components(graph)
    for i, ascc in enumerate(sccs):
        scc = order_ascc(graph, ascc.mod_ids)
        node = NodeInfo(i, scc)
        nodes.append(node)
    inv_nodes = {}  # module -> node_id
    for node in nodes:
        for mod in node.scc:
            inv_nodes[mod] = node.node_id
    for node in nodes:
        for mod in node.scc:
            state = graph[mod]
            size = 0
            if state.path:
                try:
                    size = os.path.getsize(state.path)
                except OSError:
                    pass
            node.sizes[mod] = size
            for dep in state.dependencies:
                if dep in state.priorities:
                    pri = state.priorities[dep]
                    if dep in inv_nodes:
                        dep_id = inv_nodes[dep]
                        if dep_id != node.node_id and (
                            dep_id not in node.deps or pri < node.deps[dep_id]
                        ):
                            node.deps[dep_id] = pri
    print("[" + ",\n ".join(node.dumps() for node in nodes) + "\n]", file=stdout)


def load_graph(
    sources: list[BuildSource],
    manager: BuildManager,
    old_graph: Graph | None = None,
    new_modules: list[State] | None = None,
) -> Graph:
    """Given some source files, load the full dependency graph.

    If an old_graph is passed in, it is used as the starting point and
    modified during graph loading.

    If a new_modules is passed in, any modules that are loaded are
    added to the list. This is an argument and not a return value
    so that the caller can access it even if load_graph fails.

    As this may need to parse files, this can raise CompileError in case
    there are syntax errors.
    """

    graph: Graph = old_graph if old_graph is not None else {}

    # The deque is used to implement breadth-first traversal.
    # TODO: Consider whether to go depth-first instead.  This may
    # affect the order in which we process files within import cycles.
    new = new_modules if new_modules is not None else []
    entry_points: set[str] = set()
    # Seed the graph with the initial root sources.
    for bs in sources:
        try:
            st = State.new_state(
                id=bs.module,
                path=bs.path,
                source=bs.text,
                manager=manager,
                root_source=not bs.followed,
            )
        except ModuleNotFound:
            continue
        if st.id in graph:
            manager.errors.set_file(st.xpath, st.id, manager.options)
            manager.error(
                None,
                f'Duplicate module named "{st.id}" (also at "{graph[st.id].xpath}")',
                blocker=True,
            )
            resolution_note = f"""
            See {MODULE_RESOLUTION_URL} for more info
            Common resolutions include:
                a) using `--exclude` to avoid checking one of them,
                b) adding `__init__.py` somewhere,
                c) using `--explicit-package-bases` or adjusting `MYPYPATH`
            """
            manager.note_multiline(None, resolution_note)
            manager.errors.raise_error()
        graph[st.id] = st
        new.append(st)
        entry_points.add(bs.module)
    manager.parse_all([state for state in new if state.needs_parse])

    # Note: Running this each time could be slow in the daemon. If it's a problem, we
    # can do more work to maintain this incrementally.
    seen_files = {st.abspath: st for st in graph.values() if st.path}

    # Collect dependencies.  We go breadth-first.
    # More nodes might get added to new as we go, but that's fine.
    ready = set(new)
    # Use list to make syntax error order a bit more stable.
    not_ready: list[State] = []
    for st in new:
        if st not in ready:
            # We have run out of states, parse all we have.
            assert st in not_ready
            manager.parse_all(not_ready)
            ready.update(not_ready)
            not_ready.clear()
        assert st.ancestors is not None
        # Strip out indirect dependencies.  These will be dealt with
        # when they show up as direct dependencies, and there's a
        # scenario where they hurt:
        # - Suppose A imports B and B imports C.
        # - Suppose on the next round:
        #   - C is deleted;
        #   - B is updated to remove the dependency on C;
        #   - A is unchanged.
        # - In this case A's cached *direct* dependencies are still valid
        #   (since direct dependencies reflect the imports found in the source)
        #   but A's cached *indirect* dependency on C is wrong.
        dependencies = [dep for dep in st.dependencies if st.priorities.get(dep) != PRI_INDIRECT]
        if not manager.use_fine_grained_cache():
            added = [dep for dep in st.suppressed if find_module_simple(dep, manager)]
        else:
            # During initial loading we don't care about newly added modules,
            # they will be taken care of during fine-grained update. See also
            # comment about this in `State.new_state()`.
            added = []
        for dep in st.ancestors + dependencies + st.suppressed:
            ignored = dep in st.suppressed_set and dep not in entry_points
            if ignored and dep not in added:
                manager.missing_modules[dep] = SuppressionReason.NOT_FOUND
                # TODO: for now we skip this in the daemon as a performance optimization.
                # This however creates a correctness issue, see #7777 and State.is_fresh().
                if not manager.use_fine_grained_cache() or manager.options.warn_unused_configs:
                    manager.import_options[dep] = manager.options.clone_for_module(
                        dep
                    ).dep_import_options()
            elif dep not in graph:
                try:
                    if dep in st.ancestors:
                        # TODO: Why not 'if dep not in st.dependencies' ?
                        # Ancestors don't have import context.
                        newst = State.new_state(
                            id=dep, path=None, source=None, manager=manager, ancestor_for=st
                        )
                    else:
                        newst = State.new_state(
                            id=dep,
                            path=None,
                            source=None,
                            manager=manager,
                            caller_state=st,
                            caller_line=st.dep_line_map.get(dep, 1),
                        )
                except ModuleNotFound:
                    if dep in st.dependencies_set:
                        st.suppress_dependency(dep)
                else:
                    if newst.path:
                        newst_path = newst.abspath

                        if newst_path in seen_files:
                            manager.errors.set_file(newst.xpath, newst.id, manager.options)
                            manager.error(
                                None,
                                "Source file found twice under different module names: "
                                f'"{seen_files[newst_path].id}" and "{newst.id}"',
                                blocker=True,
                            )
                            resolution_note = f"""
                            See {MODULE_RESOLUTION_URL} for more info
                            Common resolutions include:
                                a) adding `__init__.py` somewhere,
                                b) using `--explicit-package-bases` or adjusting `MYPYPATH`
                            """
                            manager.note_multiline(None, resolution_note)
                            manager.errors.raise_error()

                        seen_files[newst_path] = newst

                    assert newst.id not in graph, newst.id
                    graph[newst.id] = newst
                    new.append(newst)
                    if newst.needs_parse:
                        not_ready.append(newst)
                    else:
                        ready.add(newst)
    # There are two things we need to do after the initial load loop. One is up-suppress
    # modules that are back in graph. We need to do this after the loop to cover edge cases
    # like where a namespace package ancestor is shared by a typed and an untyped package.
    for st in graph.values():
        for dep in st.suppressed.copy():
            if dep in graph:
                st.add_dependency(dep)
                manager.missing_modules.pop(dep, None)
    # Second, in the initial loop we skip indirect dependencies, so to make indirect dependencies
    # behave more consistently with regular ones, we suppress them manually here (when needed).
    for st in graph.values():
        indirect = [dep for dep in st.dependencies if st.priorities.get(dep) == PRI_INDIRECT]
        for dep in indirect:
            if dep not in graph:
                st.suppress_dependency(dep)
    manager.plugin.set_modules(manager.modules)
    manager.errors.global_watcher = False
    return graph


def order_ascc_ex(graph: Graph, ascc: SCC) -> list[str]:
    """Apply extra heuristics on top of order_ascc().

    This should be used only for actual SCCs, not for "inner" SCCs
    we create recursively during ordering of the SCC. Currently, this
    has only some special handling for builtin SCC.
    """
    scc = order_ascc(graph, ascc.mod_ids)
    # Make the order of the SCC that includes 'builtins' and 'typing',
    # among other things, predictable. Various things may  break if
    # the order changes.
    if "builtins" in ascc.mod_ids:
        scc = sorted(scc, reverse=True)
        # If builtins is in the list, move it last.  (This is a bit of
        # a hack, but it's necessary because the builtins module is
        # part of a small cycle involving at least {builtins, abc,
        # typing}.  Of these, builtins must be processed last or else
        # some builtin objects will be incompletely processed.)
        scc.remove("builtins")
        scc.append("builtins")
    return scc


def verify_transitive_deps(ascc: SCC, graph: Graph, manager: BuildManager) -> str | None:
    """Verify all indirect dependencies of this SCC are still reachable via direct ones.

    Return first unreachable dependency id, or None.
    """
    for id in ascc.mod_ids:
        st = graph[id]
        assert st.meta is not None, "Must be called on fresh SCCs only"
        if st.trans_dep_hash == st.meta.trans_dep_hash:
            # Import graph unchanged, skip this module.
            continue
        for dep in st.dependencies:
            if st.priorities.get(dep) == PRI_INDIRECT:
                dep_scc_id = manager.scc_by_mod_id[dep].id
                if dep_scc_id == ascc.id:
                    continue
                if not manager.is_transitive_scc_dep(ascc.id, dep_scc_id):
                    return dep
    return None


def find_stale_sccs(
    sccs: list[SCC], graph: Graph, manager: BuildManager
) -> tuple[list[SCC], list[SCC]]:
    """Split a list of ready SCCs into stale and fresh.

    Fresh SCCs are those where:
    * We have valid cache files for all modules in the SCC.
    * There are no changes in dependencies (files removed from/added to the build).
    * The interface hashes of dependencies matches those recorded in the cache.
    * All indirect dependencies are still reachable via direct ones.
    The first and second conditions are verified by is_fresh().
    """
    stale_sccs = []
    fresh_sccs = []
    for ascc in sccs:
        stale_scc = {id for id in ascc.mod_ids if not graph[id].is_fresh()}
        fresh = not stale_scc

        # Verify that interfaces of dependencies still present in graph are up-to-date (fresh).
        stale_deps = set()
        for id in ascc.mod_ids:
            for dep in graph[id].dep_hashes:
                if dep in graph and graph[dep].interface_hash != graph[id].dep_hashes[dep]:
                    stale_deps.add(dep)
        fresh = fresh and not stale_deps

        # Verify the invariant that indirect dependencies are a subset of transitive direct
        # dependencies. Note: the case where indirect dependency is removed from the graph
        # completely is already handled above.
        stale_indirect = None
        if fresh:
            stale_indirect = verify_transitive_deps(ascc, graph, manager)
            if stale_indirect is not None:
                fresh = False

        if manager.logging_enabled:
            if fresh:
                fresh_msg = "fresh"
            elif stale_scc:
                fresh_msg = "inherently stale"
                if stale_scc != ascc.mod_ids:
                    fresh_msg += f" ({' '.join(sorted(stale_scc))})"
                if stale_deps:
                    fresh_msg += f" with stale deps ({' '.join(sorted(stale_deps))})"
            elif stale_deps:
                fresh_msg = f"stale due to deps ({' '.join(sorted(stale_deps))})"
            else:
                assert stale_indirect is not None
                fresh_msg = f"stale due to stale indirect dep(s): first {stale_indirect}"
            scc_str = " ".join(ascc.mod_ids)

        if fresh:
            if manager.tracing_enabled:
                manager.trace(f"Found {fresh_msg} SCC ({scc_str})")
            # If there is at most one file with errors we can skip the ordering to save time.
            mods_with_errors = [id for id in ascc.mod_ids if graph[id].error_lines]
            if len(mods_with_errors) <= 1:
                scc = mods_with_errors
            else:
                # Use exactly the same order as for stale SCCs for stability.
                scc = order_ascc_ex(graph, ascc)
            for id in scc:
                if graph[id].error_lines:
                    path = manager.errors.simplify_path(graph[id].xpath)
                    formatted = manager.errors.format_messages(
                        path, graph[id].error_lines, formatter=manager.error_formatter
                    )
                    manager.flush_errors(path, formatted, False)
            fresh_sccs.append(ascc)
        else:
            if manager.logging_enabled:
                size = len(ascc.mod_ids)
                if size == 1:
                    manager.log(f"Scheduling SCC singleton ({scc_str}) as {fresh_msg}")
                else:
                    manager.log(
                        "Scheduling SCC of size %d (%s) as %s" % (size, scc_str, fresh_msg)
                    )
            stale_sccs.append(ascc)
    return stale_sccs, fresh_sccs


def process_graph(graph: Graph, manager: BuildManager) -> None:
    """Process everything in dependency order."""
    if manager.workers:
        # Commit any cache writes from graph loading before workers try to read them.
        manager.commit()

    # Broadcast graph to workers before computing SCCs to save a bit of time.
    # TODO: check if we can optimize by sending only part of the graph needed for given SCC.
    # For example only send modules in the SCC and their dependencies.
    graph_message = GraphMessage(graph=graph, missing_modules=manager.missing_modules)
    buf = WriteBuffer()
    graph_message.write(buf)
    graph_data = buf.getvalue()
    manager.wait_ack()
    manager.broadcast(graph_data)

    sccs = sorted_components(graph)
    manager.log(
        "Found %d SCCs; largest has %d nodes" % (len(sccs), max(len(scc.mod_ids) for scc in sccs))
    )
    scc_by_id = {scc.id: scc for scc in sccs}
    manager.scc_by_id = scc_by_id
    manager.top_order = [scc.id for scc in sccs]
    for scc in sccs:
        for mod_id in scc.mod_ids:
            manager.scc_by_mod_id[mod_id] = scc

    # Broadcast SCC structure to the parallel workers, since they don't compute it.
    sccs_message = SccsDataMessage(sccs=sccs)
    buf = WriteBuffer()
    sccs_message.write(buf)
    sccs_data = buf.getvalue()
    manager.wait_ack()
    manager.broadcast(sccs_data)
    manager.wait_ack()

    manager.free_workers = set(range(manager.options.num_workers))

    # Prime the ready list with leaf SCCs (that have no dependencies).
    ready = []
    not_ready = set()
    for scc in sccs:
        if not scc.deps:
            ready.append(scc)
        else:
            not_ready.add(scc)

    still_working = False
    while ready or not_ready or still_working:
        stale, fresh = find_stale_sccs(ready, graph, manager)
        if stale:
            for scc in stale:
                for id in scc.mod_ids:
                    graph[id].mark_as_rechecked()
            manager.submit(graph, stale)
            still_working = True
        # We eagerly walk over fresh SCCs to reach as many stale SCCs as soon
        # as possible. Only when there are no fresh SCCs, we wait on scheduled stale ones.
        # This strategy, similar to a naive strategy in minesweeper game, will allow us
        # to leverage parallelism as much as possible.
        if fresh:
            done = fresh
        else:
            done, still_working, results = manager.wait_for_done(graph)
            # Expose the results of type-checking by workers. For in-process
            # type-checking this is already done and results should be empty here.
            if not manager.workers:
                assert not results
            for id, result in results.items():
                # Interface and implementation results may be mixed in the same batch
                # from different workers, process each one accordingly.
                if result.interface_hash is not None:
                    new_hash = bytes.fromhex(result.interface_hash)
                    if new_hash != graph[id].interface_hash:
                        graph[id].mark_interface_stale()
                        graph[id].interface_hash = new_hash
                else:
                    manager.flush_errors(
                        manager.errors.simplify_path(graph[id].xpath), result.error_lines, False
                    )
        ready = []
        for done_scc in done:
            for dependent in done_scc.direct_dependents:
                scc_by_id[dependent].not_ready_deps.discard(done_scc.id)
                if not scc_by_id[dependent].not_ready_deps:
                    not_ready.remove(scc_by_id[dependent])
                    ready.append(scc_by_id[dependent])
    manager.trace(f"Transitive deps cache size: {sys.getsizeof(manager.transitive_deps_cache)}")


def order_ascc(graph: Graph, ascc: AbstractSet[str], pri_max: int = PRI_INDIRECT) -> list[str]:
    """Come up with the ideal processing order within an SCC.

    Using the priorities assigned by all_imported_modules_in_file(),
    try to reduce the cycle to a DAG, by omitting arcs representing
    dependencies of lower priority.

    In the simplest case, if we have A <--> B where A has a top-level
    "import B" (medium priority) but B only has the reverse "import A"
    inside a function (low priority), we turn the cycle into a DAG by
    dropping the B --> A arc, which leaves only A --> B.

    If all arcs have the same priority, we fall back to sorting by
    reverse global order (the order in which modules were first
    encountered).

    The algorithm is recursive, as follows: when as arcs of different
    priorities are present, drop all arcs of the lowest priority,
    identify SCCs in the resulting graph, and apply the algorithm to
    each SCC thus found.  The recursion is bounded because at each
    recursion the spread in priorities is (at least) one less.

    In practice there are only a few priority levels (less than a
    dozen) and in the worst case we just carry out the same algorithm
    for finding SCCs N times.  Thus, the complexity is no worse than
    the complexity of the original SCC-finding algorithm -- see
    strongly_connected_components() below for a reference.
    """
    if len(ascc) == 1:
        return list(ascc)
    pri_spread = set()
    for id in ascc:
        state = graph[id]
        for dep in state.dependencies:
            if dep in ascc:
                pri = state.priorities.get(dep, PRI_HIGH)
                if pri < pri_max:
                    pri_spread.add(pri)
    if len(pri_spread) == 1:
        # Filtered dependencies are uniform -- order by global order.
        return sorted(ascc, key=lambda id: -graph[id].order)
    pri_max = max(pri_spread)
    sccs = sorted_components_inner(graph, ascc, pri_max)
    # The recursion is bounded by the len(pri_spread) check above.
    return [s for ss in sccs for s in order_ascc(graph, ss, pri_max)]


def process_fresh_modules(graph: Graph, modules: list[str], manager: BuildManager) -> None:
    """Process the modules in one group of modules from their cached data.

    This can be used to process an SCC of modules. This involves loading the tree (i.e.
    module symbol tables) from cache file and then fixing cross-references in the symbols.
    """
    t0 = time.time()
    for id in modules:
        graph[id].load_tree()
    t1 = time.time()
    for id in modules:
        graph[id].fix_cross_refs()
    t2 = time.time()
    manager.add_stats(process_fresh_time=t2 - t0, load_tree_time=t1 - t0)


def maybe_load_deps(graph: Graph, ascc: SCC, manager: BuildManager) -> None:
    """Load any missing fresh modules needed to process a stale SCC"""
    missing_sccs = set()
    sccs_to_find = ascc.deps.copy()
    while sccs_to_find:
        dep_scc = sccs_to_find.pop()
        if dep_scc in manager.done_sccs or dep_scc in missing_sccs:
            continue
        missing_sccs.add(dep_scc)
        sccs_to_find.update(manager.scc_by_id[dep_scc].deps)

    if missing_sccs:
        # Load missing SCCs from cache.
        # TODO: speed-up ordering if this causes problems for large builds.
        fresh_sccs_to_load = [
            manager.scc_by_id[sid] for sid in manager.top_order if sid in missing_sccs
        ]

        if manager.parallel_worker:
            # Update cache metas as well, cache data is loaded below
            # in process_fresh_modules().
            for prev_scc in fresh_sccs_to_load:
                for mod_id in prev_scc.mod_ids:
                    graph[mod_id].reload_meta()

        manager.log(f"Processing {len(fresh_sccs_to_load)} fresh SCCs")
        if (
            not manager.options.test_env
            and platform.python_implementation() == "CPython"
            and manager.gc_freeze_cycles < MAX_GC_FREEZE_CYCLES
        ):
            # When deserializing cache we create huge amount of new objects, so even
            # with our generous GC thresholds, GC is still doing a lot of pointless
            # work searching for garbage. So, we temporarily disable it when
            # processing fresh SCCs, and then move all the new objects to the oldest
            # generation with the freeze()/unfreeze() trick below. This is arguably
            # a hack, but it gives huge performance wins for large third-party
            # libraries, like torch.
            gc.disable()
        for prev_scc in fresh_sccs_to_load:
            manager.done_sccs.add(prev_scc.id)
            process_fresh_modules(graph, sorted(prev_scc.mod_ids), manager)
        if (
            not manager.options.test_env
            and platform.python_implementation() == "CPython"
            and manager.gc_freeze_cycles < MAX_GC_FREEZE_CYCLES
        ):
            manager.gc_freeze_cycles += 1
            gc.freeze()
            gc.unfreeze()
            gc.enable()


def process_stale_scc(graph: Graph, ascc: SCC, manager: BuildManager) -> None:
    """Process the modules in one SCC from source code."""
    # First verify if all transitive dependencies are loaded in the current process.
    t0 = time.time()
    maybe_load_deps(graph, ascc, manager)
    t1 = time.time()
    # Process the SCC in stable order.
    scc = order_ascc_ex(graph, ascc)

    t2 = time.time()
    stale = scc
    # Parse before verify_dependencies so that inline config comments
    # (e.g. "# mypy: disable-error-code") are applied to options.
    manager.parse_all([graph[id] for id in stale], post_parse=False)
    for id in stale:
        # Re-generate import errors in case this module was loaded from the cache.
        if graph[id].meta:
            graph[id].verify_dependencies(suppressed_only=True)
    if "typing" in scc:
        # For historical reasons we need to manually add typing aliases
        # for built-in generic collections, see docstring of
        # SemanticAnalyzerPass2.add_builtin_aliases for details.
        typing_mod = graph["typing"].tree
        assert typing_mod, "The typing module was not parsed"
    mypy.semanal_main.semantic_analysis_for_scc(graph, scc, manager.errors)

    t3 = time.time()
    # Track what modules aren't yet done, so we can finish them as soon
    # as possible, saving memory.
    unfinished_modules = set(stale)
    for id in stale:
        graph[id].type_check_first_pass()
        if not graph[id].type_checker().deferred_nodes:
            unfinished_modules.discard(id)
            graph[id].detect_possibly_undefined_vars()
            graph[id].finish_passes()

    while unfinished_modules:
        for id in stale:
            if id not in unfinished_modules:
                continue
            if not graph[id].type_check_second_pass():
                unfinished_modules.discard(id)
                graph[id].detect_possibly_undefined_vars()
                graph[id].finish_passes()
    for id in stale:
        graph[id].generate_unused_ignore_notes()
        graph[id].generate_ignore_without_code_notes()

    t4 = time.time()
    # Flush errors, and write cache in two phases: first data files, then meta files.
    # The two-phase structure is needed because meta.dep_hashes references interface_hash
    # values from other modules in the SCC, which are updated by write_cache().
    meta_tuples = {}
    errors_by_id = {}
    for id in stale:
        if graph[id].xpath not in manager.errors.ignored_files:
            errors = manager.errors.file_messages(graph[id].xpath)
            formatted = manager.errors.format_messages(
                graph[id].xpath, errors, formatter=manager.error_formatter
            )
            manager.flush_errors(manager.errors.simplify_path(graph[id].xpath), formatted, False)
            errors_by_id[id] = errors
        meta_tuple = graph[id].write_cache()
        meta_tuples[id] = meta_tuple
        # Commit data file write immediately to avoid holding shard locks across modules.
        if meta_tuple is not None:
            manager.commit_module(meta_tuple[1])
    for id in stale:
        meta_tuple = meta_tuples[id]
        if meta_tuple is None:
            continue
        meta, meta_file = meta_tuple
        state = graph[id]
        # Indirect dependencies are stored as part of CacheMetaEx below.
        meta.dep_hashes = [
            graph[dep].interface_hash
            for dep in graph[id].dependencies
            if state.priorities.get(dep) != PRI_INDIRECT
        ]
        write_cache_meta(meta, manager, meta_file)
        indirect = [dep for dep in state.dependencies if state.priorities.get(dep) == PRI_INDIRECT]
        meta_ex = CacheMetaEx(
            dependencies=indirect,
            suppressed=[
                dep for dep in state.suppressed if state.priorities.get(dep) == PRI_INDIRECT
            ],
            dep_hashes=[graph[dep].interface_hash for dep in indirect],
            error_lines=errors_by_id.get(id, []),
        )
        write_cache_meta_ex(meta_file, meta_ex, manager)
        manager.commit_module(meta_file)
    manager.done_sccs.add(ascc.id)
    manager.add_stats(
        load_missing_time=t1 - t0,
        order_scc_time=t2 - t1,
        semanal_time=t3 - t2,
        type_check_time=t4 - t3,
        flush_and_cache_time=time.time() - t4,
    )


def process_stale_scc_interface(
    graph: Graph, ascc: SCC, manager: BuildManager, from_cache: set[str]
) -> list[tuple[str, ModuleResult, str]]:
    """Process the modules' interfaces in one SCC from source code."""
    # First verify if all transitive dependencies are loaded in the current process.
    t0 = time.time()
    maybe_load_deps(graph, ascc, manager)
    t1 = time.time()
    # Process the SCC in stable order.
    scc = order_ascc_ex(graph, ascc)

    t2 = time.time()
    stale = scc
    for id in stale:
        # Re-generate import errors in case this module was loaded from the cache.
        # Deserialized states all have meta=None, so the caller should specify
        # explicitly which of them are from cache.
        if id in from_cache:
            graph[id].verify_dependencies(suppressed_only=True)
    mypy.semanal_main.semantic_analysis_for_scc(graph, scc, manager.errors)

    t3 = time.time()
    # Track what modules aren't yet done, so we can finish them as soon
    # as possible, saving memory.
    unfinished_modules = set(stale)
    for id in stale:
        graph[id].type_check_first_pass(recurse_into_functions=False)
        if not graph[id].type_checker().deferred_nodes:
            unfinished_modules.discard(id)
    while unfinished_modules:
        for id in stale:
            if id not in unfinished_modules:
                continue
            if not graph[id].type_check_second_pass(recurse_into_functions=False):
                unfinished_modules.discard(id)

    t4 = time.time()
    scc_result = []
    meta_tuples = {}
    for id in stale:
        meta_tuple = graph[id].write_cache()
        meta_tuples[id] = meta_tuple
        # Commit data file write immediately to avoid holding shard locks across modules.
        if meta_tuple is not None:
            manager.commit_module(meta_tuple[1])
    for id in stale:
        meta_tuple = meta_tuples[id]
        if meta_tuple is None:
            continue
        meta, meta_file = meta_tuple
        state = graph[id]
        meta.dep_hashes = [
            graph[dep].interface_hash
            for dep in state.dependencies
            if state.priorities.get(dep) != PRI_INDIRECT
        ]
        write_cache_meta(meta, manager, meta_file)
        manager.commit_module(meta_file)
        scc_result.append((id, ModuleResult(graph[id].interface_hash.hex(), []), meta_file))
    manager.done_sccs.add(ascc.id)
    manager.add_stats(
        load_missing_time=t1 - t0,
        order_scc_time=t2 - t1,
        semanal_time=t3 - t2,
        type_check_time_interface=t4 - t3,
        flush_and_cache_time=time.time() - t4,
    )
    return scc_result


def process_stale_scc_implementation(
    graph: Graph, stale: list[str], manager: BuildManager, meta_files: list[str]
) -> dict[str, ModuleResult]:
    """Process implementations (top-level function/method bodies) in an SCC."""
    t0 = time.time()
    unfinished_modules = set(stale)
    for id in stale:
        checker = graph[id].type_checker()
        # Optimization: if this is a 3rd party library, or we ignore errors
        # otherwise in this module, skip the implementations altogether.
        if checker.can_skip_diagnostics and not checker.options.preserve_asts:
            unfinished_modules.discard(id)
            graph[id].finish_passes()
            continue
        # We need to reset deferral count after possibly deferring any methods that
        # are considered part of the top-level (because they define/infer variables).
        checker.pass_num = 0
        checker.deferred_nodes.clear()
        tree = graph[id].tree
        assert tree is not None
        todo = []
        # Passing impl_only will select only "leaf" nodes (not the TypeInfos).
        for _, node, info in tree.local_definitions(impl_only=True):
            assert isinstance(node.node, (FuncDef, OverloadedFuncDef, Decorator))
            todo.append(DeferredNode(node.node, info))
        graph[id].type_check_second_pass(todo=todo, impl_only=True)
        if not checker.deferred_nodes:
            unfinished_modules.discard(id)
            graph[id].detect_possibly_undefined_vars()
            graph[id].finish_passes()
    while unfinished_modules:
        for id in stale:
            if id not in unfinished_modules:
                continue
            if not graph[id].type_check_second_pass(impl_only=True):
                unfinished_modules.discard(id)
                graph[id].detect_possibly_undefined_vars()
                graph[id].finish_passes()

    for id in stale:
        graph[id].generate_unused_ignore_notes()
        graph[id].generate_ignore_without_code_notes()

    scc_result = {}
    for id, meta_file in zip(stale, meta_files):
        state = graph[id]
        indirect = [dep for dep in state.dependencies if state.priorities.get(dep) == PRI_INDIRECT]
        meta_ex = CacheMetaEx(
            dependencies=indirect,
            suppressed=[
                dep for dep in state.suppressed if state.priorities.get(dep) == PRI_INDIRECT
            ],
            dep_hashes=[graph[dep].interface_hash for dep in indirect],
            error_lines=[],
        )
        if graph[id].xpath not in manager.errors.ignored_files:
            errors = manager.errors.file_messages(graph[id].xpath)
            formatted = manager.errors.format_messages(
                graph[id].xpath, errors, formatter=manager.error_formatter
            )
            meta_ex.error_lines = errors
            write_cache_meta_ex(meta_file, meta_ex, manager)
            scc_result[id] = ModuleResult(None, formatted)
        else:
            # If there are no errors, only write the cache, don't send anything back
            # to the caller (as a micro-optimization).
            write_cache_meta_ex(meta_file, meta_ex, manager)
        manager.commit_module(meta_file)

    manager.add_stats(type_check_time_implementation=time.time() - t0)
    return scc_result


def prepare_sccs_full(
    raw_sccs: Iterator[set[str]], edges: dict[str, list[str]]
) -> dict[SCC, set[SCC]]:
    """Turn raw SCC sets into SCC objects and build dependency graph for SCCs."""
    sccs = [SCC(raw_scc) for raw_scc in raw_sccs]
    scc_map = {}
    for scc in sccs:
        for id in scc.mod_ids:
            scc_map[id] = scc
    scc_deps_map: dict[SCC, set[SCC]] = {}
    for scc in sccs:
        for id in scc.mod_ids:
            scc_deps_map.setdefault(scc, set()).update(scc_map[dep] for dep in edges[id])
    for scc in sccs:
        # Remove trivial dependency on itself.
        scc_deps_map[scc].discard(scc)
        for dep_scc in scc_deps_map[scc]:
            scc.deps.add(dep_scc.id)
            scc.not_ready_deps.add(dep_scc.id)
    return scc_deps_map


def sorted_components(graph: Graph) -> list[SCC]:
    """Return the graph's SCCs, topologically sorted by dependencies.

    The sort order is from leaves (nodes without dependencies) to
    roots (nodes on which no other nodes depend).
    """
    # Compute SCCs.
    vertices = set(graph)
    edges = {id: deps_filtered(graph, vertices, id, PRI_INDIRECT) for id in vertices}
    scc_dep_map = prepare_sccs_full(strongly_connected_components(vertices, edges), edges)
    # Topsort.
    res = []
    for ready in topsort(scc_dep_map):
        # Sort the sets in ready by reversed smallest State.order.  Examples:
        #
        # - If ready is [{x}, {y}], x.order == 1, y.order == 2, we get
        #   [{y}, {x}].
        #
        # - If ready is [{a, b}, {c, d}], a.order == 1, b.order == 3,
        #   c.order == 2, d.order == 4, the sort keys become [1, 2]
        #   and the result is [{c, d}, {a, b}].
        sorted_ready = sorted(ready, key=lambda scc: -min(graph[id].order for id in scc.mod_ids))
        for scc in sorted_ready:
            scc.size_hint = sum(graph[mid].size_hint for mid in scc.mod_ids)
            for dep in scc_dep_map[scc]:
                dep.direct_dependents.append(scc.id)
            # We compute dependencies hash here since we know no direct
            # dependencies will be added or suppressed after this point.
            trans_dep_hash = transitive_dep_hash(scc, graph)
            for id in scc.mod_ids:
                graph[id].trans_dep_hash = trans_dep_hash
        res.extend(sorted_ready)
    return res


def sorted_components_inner(
    graph: Graph, vertices: AbstractSet[str], pri_max: int
) -> list[AbstractSet[str]]:
    """Simplified version of sorted_components() to work with sub-graphs.

    This doesn't create SCC objects, and operates with raw sets. This function
    also allows filtering dependencies to take into account when building SCCs.
    This is used for heuristic ordering of modules within actual SCCs.
    """
    edges = {id: deps_filtered(graph, vertices, id, pri_max) for id in vertices}
    sccs = list(strongly_connected_components(vertices, edges))
    res = []
    for ready in topsort(prepare_sccs(sccs, edges)):
        res.extend(sorted(ready, key=lambda scc: -min(graph[id].order for id in scc)))
    return res


def deps_filtered(graph: Graph, vertices: AbstractSet[str], id: str, pri_max: int) -> list[str]:
    """Filter dependencies for id with pri < pri_max."""
    if id not in vertices:
        return []
    state = graph[id]
    return [
        dep
        for dep in state.dependencies
        if dep in vertices and state.priorities.get(dep, PRI_HIGH) < pri_max
    ]


def transitive_dep_hash(scc: SCC, graph: Graph) -> bytes:
    """Compute stable snapshot of transitive import structure for given SCC."""
    all_direct_deps = sorted(
        {
            dep
            for id in scc.mod_ids
            for dep in graph[id].dependencies
            if graph[id].priorities.get(dep) != PRI_INDIRECT
        }
    )
    buf = WriteBuffer()
    for dep_id in all_direct_deps:
        write_str_bare(buf, dep_id)
        if dep_id not in scc.mod_ids:
            write_bytes_bare(buf, graph[dep_id].trans_dep_hash)
    return hash_digest_bytes(buf.getvalue())


def missing_stubs_file(cache_dir: str) -> str:
    return os_path_join(cache_dir, "missing_stubs")


def record_missing_stub_packages(cache_dir: str, missing_stub_packages: set[str]) -> None:
    """Write a file containing missing stub packages.

    This allows a subsequent "mypy --install-types" run (without other arguments)
    to install missing stub packages.
    """
    fnam = missing_stubs_file(cache_dir)
    if missing_stub_packages:
        with open(fnam, "w") as f:
            for pkg in sorted(missing_stub_packages):
                f.write(f"{pkg}\n")
    else:
        if os.path.isfile(fnam):
            os.remove(fnam)


def is_silent_import_module(manager: BuildManager, path: str) -> bool:
    if manager.options.no_silence_site_packages:
        return False
    # Silence errors in site-package dirs and typeshed
    if any(is_sub_path_normabs(path, dir) for dir in manager.search_paths.package_path):
        return True
    return any(is_sub_path_normabs(path, dir) for dir in manager.search_paths.typeshed_path)


def write_undocumented_ref_info(
    state: State, metastore: MetadataStore, options: Options, type_map: dict[Expression, Type]
) -> None:
    # This exports some dependency information in a rather ad-hoc fashion, which
    # can be helpful for some tools. This is all highly experimental and could be
    # removed at any time.

    from mypy.refinfo import get_undocumented_ref_info_json

    if not state.tree:
        # We need a full AST for this.
        return

    _, data_file, _ = get_cache_names(state.id, state.xpath, options)
    ref_info_file = ".".join(data_file.split(".")[:-2]) + ".refs.json"
    assert not ref_info_file.startswith(".")

    deps_json = get_undocumented_ref_info_json(state.tree, type_map)
    metastore.write(ref_info_file, json_dumps(deps_json))


# The IPC message classes and tags for communication with build workers are
# in this file to avoid import cycles.
# Note that we use a more compact fixed serialization format than in cache.py.
# This is because the messages don't need to read by a generic tool, nor there
# is any need for backwards compatibility. We still reuse some elements from
# cache.py for convenience, and also some conventions (like using bare ints
# to specify object size).
# Note that we can use tags overlapping with cache.py, since they should never
# appear on the same context.
ACK_MESSAGE: Final[Tag] = 101
SCC_REQUEST_MESSAGE: Final[Tag] = 102
SCC_RESPONSE_MESSAGE: Final[Tag] = 103
SOURCES_DATA_MESSAGE: Final[Tag] = 104
SCCS_DATA_MESSAGE: Final[Tag] = 105
GRAPH_MESSAGE: Final[Tag] = 106


class AckMessage(IPCMessage):
    """An empty message used primarily for synchronization."""

    @classmethod
    def read(cls, buf: ReadBuffer) -> AckMessage:
        return AckMessage()

    def write(self, buf: WriteBuffer) -> None:
        write_tag(buf, ACK_MESSAGE)


class SccRequestMessage(IPCMessage):
    """
    A message representing a request to type check a batch of SCCs.

    If scc_ids is empty, then it means that the coordinator requested a shutdown.
    """

    def __init__(
        self,
        *,
        scc_ids: list[int],
        import_errors: dict[str, list[ErrorInfo]],
        mod_data: dict[str, tuple[bytes, FileRawData | None]],
    ) -> None:
        self.scc_ids = scc_ids
        self.import_errors = import_errors
        self.mod_data = mod_data

    @classmethod
    def read(cls, buf: ReadBuffer) -> SccRequestMessage:
        return SccRequestMessage(
            scc_ids=read_int_list(buf),
            import_errors={
                read_str(buf): [ErrorInfo.read(buf) for _ in range(read_int_bare(buf))]
                for _ in range(read_int_bare(buf))
            },
            mod_data={
                read_str_bare(buf): (
                    read_bytes(buf),
                    FileRawData.read(buf) if read_bool(buf) else None,
                )
                for _ in range(read_int_bare(buf))
            },
        )

    def write(self, buf: WriteBuffer) -> None:
        write_tag(buf, SCC_REQUEST_MESSAGE)
        write_int_list(buf, self.scc_ids)
        write_int_bare(buf, len(self.import_errors))
        for path, errors in self.import_errors.items():
            write_str(buf, path)
            write_int_bare(buf, len(errors))
            for error in errors:
                error.write(buf)
        write_int_bare(buf, len(self.mod_data))
        for mod, (suppressed_deps_opts, raw_data) in self.mod_data.items():
            write_str_bare(buf, mod)
            write_bytes(buf, suppressed_deps_opts)
            if raw_data is None:
                write_bool(buf, False)
            else:
                write_bool(buf, True)
                raw_data.write(buf)


class ModuleResult:
    """A simple class representing the result of type-checking phase for a single module.

    Non-None interface hash signifies this is a result of checking the interface
    of the module, otherwise this is a result of checking the implementation (which
    includes errors encountered during both phases).
    """

    def __init__(self, interface_hash: str | None, error_lines: list[str]) -> None:
        self.interface_hash = interface_hash
        self.error_lines = error_lines

    @classmethod
    def read(cls, buf: ReadBuffer) -> ModuleResult:
        return ModuleResult(read_str_opt(buf), read_str_list(buf))

    def write(self, buf: WriteBuffer) -> None:
        write_str_opt(buf, self.interface_hash)
        write_str_list(buf, self.error_lines)


class SccResponseMessage(IPCMessage):
    """
    A message representing a result of type checking a batch of SCCs.

    Only one of `result` or `blocker` can be non-None. The latter means there was
    a blocking error while type checking the SCCs. The `is_interface` flag indicates
    whether this is a result for interface or implementation phase of type-checking.
    """

    def __init__(
        self,
        *,
        scc_ids: list[int],
        is_interface: bool,
        result: dict[str, ModuleResult] | None = None,
        blocker: CompileError | None = None,
    ) -> None:
        if result is not None:
            assert blocker is None
        if blocker is not None:
            assert result is None
        self.scc_ids = scc_ids
        self.is_interface = is_interface
        self.result = result
        self.blocker = blocker

    @classmethod
    def read(cls, buf: ReadBuffer) -> SccResponseMessage:
        scc_ids = read_int_list(buf)
        is_interface = read_bool(buf)
        tag = read_tag(buf)
        if tag == LITERAL_NONE:
            return SccResponseMessage(
                scc_ids=scc_ids,
                is_interface=is_interface,
                blocker=CompileError(read_str_list(buf), read_bool(buf), read_str_opt(buf)),
            )
        else:
            assert tag == DICT_STR_GEN
            return SccResponseMessage(
                scc_ids=scc_ids,
                is_interface=is_interface,
                result={
                    read_str_bare(buf): ModuleResult.read(buf) for _ in range(read_int_bare(buf))
                },
            )

    def write(self, buf: WriteBuffer) -> None:
        write_tag(buf, SCC_RESPONSE_MESSAGE)
        write_int_list(buf, self.scc_ids)
        write_bool(buf, self.is_interface)
        if self.result is None:
            assert self.blocker is not None
            write_tag(buf, LITERAL_NONE)
            write_str_list(buf, self.blocker.messages)
            write_bool(buf, self.blocker.use_stdout)
            write_str_opt(buf, self.blocker.module_with_blocker)
        else:
            write_tag(buf, DICT_STR_GEN)
            write_int_bare(buf, len(self.result))
            for mod_id in sorted(self.result):
                write_str_bare(buf, mod_id)
                self.result[mod_id].write(buf)


class SourcesDataMessage(IPCMessage):
    """A message wrapping a list of build sources."""

    def __init__(self, *, sources: list[BuildSource]) -> None:
        self.sources = sources

    @classmethod
    def read(cls, buf: ReadBuffer) -> SourcesDataMessage:
        sources = [
            BuildSource(
                read_str_opt(buf),
                read_str_opt(buf),
                read_str_opt(buf),
                read_str_opt(buf),
                read_bool(buf),
            )
            for _ in range(read_int_bare(buf))
        ]
        return SourcesDataMessage(sources=sources)

    def write(self, buf: WriteBuffer) -> None:
        write_tag(buf, SOURCES_DATA_MESSAGE)
        write_int_bare(buf, len(self.sources))
        for bs in self.sources:
            write_str_opt(buf, bs.path)
            write_str_opt(buf, bs.module)
            write_str_opt(buf, bs.text)
            write_str_opt(buf, bs.base_dir)
            write_bool(buf, bs.followed)


class SccsDataMessage(IPCMessage):
    """A message wrapping the SCC structure computed by the coordinator."""

    def __init__(self, *, sccs: list[SCC]) -> None:
        self.sccs = sccs

    @classmethod
    def read(cls, buf: ReadBuffer) -> SccsDataMessage:
        sccs = [
            SCC(set(read_str_list(buf)), read_int(buf), read_int_list(buf))
            for _ in range(read_int_bare(buf))
        ]
        return SccsDataMessage(sccs=sccs)

    def write(self, buf: WriteBuffer) -> None:
        write_tag(buf, SCCS_DATA_MESSAGE)
        write_int_bare(buf, len(self.sccs))
        for scc in self.sccs:
            write_str_list(buf, sorted(scc.mod_ids))
            write_int(buf, scc.id)
            write_int_list(buf, sorted(scc.deps))


class GraphMessage(IPCMessage):
    """A message wrapping the build graph computed by the coordinator."""

    def __init__(self, *, graph: Graph, missing_modules: dict[str, int]) -> None:
        self.graph = graph
        self.missing_modules = missing_modules
        # Send this data separately as it will be lost during state serialization.
        self.from_cache = {mod_id for mod_id in graph if graph[mod_id].meta}

    @classmethod
    def read(cls, buf: ReadBuffer, manager: BuildManager | None = None) -> GraphMessage:
        assert manager is not None
        graph = {read_str_bare(buf): State.read(buf, manager) for _ in range(read_int_bare(buf))}
        missing_modules = {read_str_bare(buf): read_int(buf) for _ in range(read_int_bare(buf))}
        message = GraphMessage(graph=graph, missing_modules=missing_modules)
        message.from_cache = {read_str_bare(buf) for _ in range(read_int_bare(buf))}
        return message

    def write(self, buf: WriteBuffer) -> None:
        write_tag(buf, GRAPH_MESSAGE)
        write_int_bare(buf, len(self.graph))
        for mod_id, state in self.graph.items():
            write_str_bare(buf, mod_id)
            state.write(buf)
        write_int_bare(buf, len(self.missing_modules))
        for module, reason in self.missing_modules.items():
            write_str_bare(buf, module)
            write_int(buf, reason)
        write_int_bare(buf, len(self.from_cache))
        for module in self.from_cache:
            write_str_bare(buf, module)
