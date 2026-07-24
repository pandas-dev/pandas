"""
Mypy parallel build worker.

The protocol of communication with the coordinator is as following:
* Read (pickled) build options from command line.
* Populate status file with pid and socket address.
* Receive build sources from coordinator.
* Initialize build manager with the sources, and send ack to coordinator.
* Receive build graph from coordinator, and ack it.
* Receive SCC structure from coordinator, and ack it.
* In a loop:
  - Receive an SCC id from coordinator, and start processing it.
  - SCC is processed in two phases: interface and implementation, send a response after each.
  - When prompted by coordinator (with a scc_ids=[] message), cleanup and shutdown.
"""

from __future__ import annotations

import argparse
import gc
import json
import os
import platform
import sys
import time
from typing import Final

from librt.internal import ReadBuffer, read_tag

from mypy import util
from mypy.build import (
    GRAPH_MESSAGE,
    SCC_REQUEST_MESSAGE,
    SCCS_DATA_MESSAGE,
    SOURCES_DATA_MESSAGE,
    AckMessage,
    BuildManager,
    Graph,
    GraphMessage,
    ModuleResult,
    SccRequestMessage,
    SccResponseMessage,
    SccsDataMessage,
    SourcesDataMessage,
    load_plugins,
    process_stale_scc_implementation,
    process_stale_scc_interface,
)
from mypy.cache import Tag, read_int_list, read_json
from mypy.defaults import RECURSION_LIMIT, WORKER_CONNECTION_TIMEOUT, WORKER_IDLE_TIMEOUT
from mypy.error_formatter import OUTPUT_CHOICES
from mypy.errors import CompileError, ErrorInfo, Errors, report_internal_error
from mypy.fscache import FileSystemCache
from mypy.ipc import IPCException, IPCServer, ready_to_read, receive, send
from mypy.modulefinder import BuildSource, BuildSourceSet, compute_search_paths
from mypy.nodes import FileRawData
from mypy.options import Options
from mypy.util import read_py_file
from mypy.version import __version__

parser = argparse.ArgumentParser(prog="mypy_worker", description="Mypy build worker")
parser.add_argument("--status-file", help="status file to communicate worker details")
parser.add_argument("--options-data", help="file with serialized mypy options")

CONNECTION_NAME = "build_worker"


class ServerContext:
    def __init__(
        self,
        options: Options,
        disable_error_code: list[str],
        enable_error_code: list[str],
        errors: Errors,
        fscache: FileSystemCache,
    ) -> None:
        self.options: Final = options
        self.disable_error_code: Final = disable_error_code
        self.enable_error_code: Final = enable_error_code
        self.errors: Final = errors
        self.fscache: Final = fscache


def main(argv: list[str]) -> None:
    # Set recursion limit and GC thresholds consistent with mypy/main.py
    sys.setrecursionlimit(RECURSION_LIMIT)
    if platform.python_implementation() == "CPython":
        gc.set_threshold(200 * 1000, 30, 30)

    args = parser.parse_args(argv)

    # This mimics how daemon receives the options. Note we need to postpone
    # processing error codes after plugins are loaded, because plugins can add
    # custom error codes.
    with open(args.options_data, "rb") as f:
        buf = ReadBuffer(f.read())
    options_dict = read_json(buf)
    disable_error_code = options_dict.pop("disable_error_code", [])
    enable_error_code = options_dict.pop("enable_error_code", [])
    options = Options().apply_changes(options_dict)

    status_file = args.status_file
    server = IPCServer(CONNECTION_NAME, WORKER_CONNECTION_TIMEOUT)

    try:
        with open(status_file, "w") as f:
            json.dump({"pid": os.getpid(), "connection_name": server.connection_name}, f)
            f.write("\n")
    except Exception as exc:
        print(f"Error writing status file {status_file}:", exc)
        raise

    fscache = FileSystemCache()
    fscache.set_package_root(options.package_root)
    cached_read = fscache.read
    error_formatter = None if options.output is None else OUTPUT_CHOICES.get(options.output)
    errors = Errors(
        options,
        read_source=lambda path: read_py_file(path, cached_read),
        error_formatter=error_formatter,
    )

    ctx = ServerContext(options, disable_error_code, enable_error_code, errors, fscache)
    try:
        with server:
            serve(server, ctx)
    except (OSError, IPCException) as exc:
        if options.verbosity >= 1:
            print("Error communicating with coordinator:", exc)
    except Exception as exc:
        report_internal_error(exc, errors.file, 0, errors, options)
    finally:
        server.cleanup()

    if options.fast_exit:
        # Exit fast if allowed, since coordinator is waiting on us.
        util.hard_exit(0)


def should_shutdown(buf: ReadBuffer, expected_tag: Tag) -> bool:
    """Check if the message is a shutdown request."""
    tag = read_tag(buf)
    if tag == SCC_REQUEST_MESSAGE:
        assert not read_int_list(buf)
        return True
    assert tag == expected_tag, f"Unexpected tag: {tag}"
    return False


def serve(server: IPCServer, ctx: ServerContext) -> None:
    """Main server loop of the worker.

    Receive initial state from the coordinator, then process each
    SCC checking request and reply to client (coordinator). See module
    docstring for more details on the protocol.
    """
    buf = receive(server)
    if should_shutdown(buf, SOURCES_DATA_MESSAGE):
        return
    sources = SourcesDataMessage.read(buf).sources
    manager = setup_worker_manager(sources, ctx)
    if manager is None:
        return

    # Notify coordinator we are done with setup.
    send(server, AckMessage())
    buf = receive(server)
    if should_shutdown(buf, GRAPH_MESSAGE):
        return

    # Disable GC before loading graph and SCC structure, these create a bunch
    # of small objects that will stay around until the end of the build.
    if platform.python_implementation() == "CPython":
        gc.disable()

    graph_data = GraphMessage.read(buf, manager)
    # Update some manager data in-place as it has been passed to semantic analyzer.
    manager.missing_modules |= graph_data.missing_modules
    graph = graph_data.graph
    for id in graph:
        manager.import_map[id] = graph[id].dependencies_set
    # Link modules dicts, so that plugins will get access to ASTs as we parse them.
    manager.plugin.set_modules(manager.modules)

    # Notify coordinator we are ready to receive computed graph SCC structure.
    send(server, AckMessage())
    buf = receive(server)
    if should_shutdown(buf, SCCS_DATA_MESSAGE):
        return
    sccs = SccsDataMessage.read(buf).sccs
    manager.scc_by_id = {scc.id: scc for scc in sccs}
    manager.top_order = [scc.id for scc in sccs]

    if platform.python_implementation() == "CPython":
        gc.freeze()
        gc.enable()

    # Notify coordinator we are ready to start processing SCCs.
    send(server, AckMessage())
    while True:
        t0 = time.time()
        ready_to_read([server], WORKER_IDLE_TIMEOUT)
        t1 = time.time()
        buf = receive(server)
        assert read_tag(buf) == SCC_REQUEST_MESSAGE
        scc_message = SccRequestMessage.read(buf)
        manager.add_stats(scc_wait_time=t1 - t0, scc_receive_time=time.time() - t1)
        scc_ids = scc_message.scc_ids
        if not scc_ids:
            # This indicates a shutdown request. Add GC stats before exiting.
            gc_stats = gc.get_stats()
            manager.add_stats(
                gc_collections_gen0=gc_stats[0]["collections"],
                gc_collections_gen1=gc_stats[1]["collections"],
                gc_collections_gen2=gc_stats[2]["collections"],
            )
            manager.dump_stats()
            break
        sccs = [manager.scc_by_id[scc_id] for scc_id in scc_ids]
        mod_ids: list[str] = []
        for scc in sccs:
            mod_ids.extend(scc.mod_ids)
        t0 = time.time()
        try:
            load_states(mod_ids, graph, manager, scc_message.import_errors, scc_message.mod_data)
            results = []
            for scc in sccs:
                scc_result = process_stale_scc_interface(
                    graph, scc, manager, from_cache=graph_data.from_cache
                )
                results.extend(scc_result)
                # We must commit after each SCC, otherwise we break --sqlite-cache.
                manager.commit()
        except CompileError as blocker:
            message = SccResponseMessage(scc_ids=scc_ids, is_interface=True, blocker=blocker)
            timed_send(manager, server, message)
        else:
            mod_results = {}
            stale = []
            meta_files = []
            for id, mod_result, meta_file in results:
                stale.append(id)
                mod_results[id] = mod_result
                meta_files.append(meta_file)
            message = SccResponseMessage(scc_ids=scc_ids, is_interface=True, result=mod_results)
            timed_send(manager, server, message)
            try:
                # Process implementations one by one, so that we can free memory a bit earlier.
                result: dict[str, ModuleResult] = {}
                for id, meta_file in zip(stale, meta_files):
                    result |= process_stale_scc_implementation(graph, [id], manager, [meta_file])
                # Both phases write cache, so we should commit here as well.
                manager.commit()
            except CompileError as blocker:
                message = SccResponseMessage(scc_ids=scc_ids, is_interface=False, blocker=blocker)
            else:
                message = SccResponseMessage(scc_ids=scc_ids, is_interface=False, result=result)
            timed_send(manager, server, message)
        manager.add_stats(total_process_stale_time=time.time() - t0, stale_sccs_processed=1)


def timed_send(manager: BuildManager, server: IPCServer, message: SccResponseMessage) -> None:
    t0 = time.time()
    send(server, message)
    manager.add_stats(scc_send_time=time.time() - t0)


def load_states(
    mod_ids: list[str],
    graph: Graph,
    manager: BuildManager,
    import_errors: dict[str, list[ErrorInfo]],
    mod_data: dict[str, tuple[bytes, FileRawData | None]],
) -> None:
    """Re-create full state of an SCC as it would have been in coordinator."""
    if platform.python_implementation() == "CPython":
        # Run full collection after previous SCC batch, everything that survives
        # will be put into permanent generation below, since we don't free anything
        # after SCC processing is done.
        gc.collect()
        gc.disable()
    needs_parse = []
    for id in mod_ids:
        state = graph[id]
        # Re-clone options since we don't send them, it is usually faster than deserializing.
        state.options = state.options.clone_for_module(state.id)
        suppressed_deps_opts, raw_data = mod_data[id]
        if raw_data is not None:
            state.parse_file(raw_data=raw_data)
        else:
            needs_parse.append(state)
        # Set data that is needed to be written to cache meta.
        state.known_suppressed_deps_opts = suppressed_deps_opts
    # Perform actual parsing in parallel (but we don't need to compute dependencies).
    if needs_parse:
        manager.parse_all(needs_parse, post_parse=False)
    for id in mod_ids:
        state = graph[id]
        assert state.tree is not None
        import_lines = {imp.line for imp in state.tree.imports}
        state.imports_ignored = {
            line: codes for line, codes in state.tree.ignored_lines.items() if line in import_lines
        }
        # Replay original errors encountered during graph loading in coordinator.
        if id in import_errors:
            manager.errors.set_file(state.xpath, id, state.options)
            for err_info in import_errors[id]:
                manager.errors.add_error_info(err_info)
    if platform.python_implementation() == "CPython":
        gc.freeze()
        gc.enable()


def setup_worker_manager(sources: list[BuildSource], ctx: ServerContext) -> BuildManager | None:
    data_dir = os.path.dirname(os.path.dirname(__file__))
    # This is used for testing only now.
    alt_lib_path = os.environ.get("MYPY_ALT_LIB_PATH")
    search_paths = compute_search_paths(sources, ctx.options, data_dir, alt_lib_path)

    source_set = BuildSourceSet(sources)
    try:
        plugin, snapshot = load_plugins(ctx.options, ctx.errors, sys.stdout, [])
    except CompileError:
        # CompileError while importing plugins will be reported by the coordinator.
        return None

    # Process the rest of the options when plugins are loaded.
    options = ctx.options
    options.disable_error_code = ctx.disable_error_code
    options.enable_error_code = ctx.enable_error_code
    options.process_error_codes(error_callback=lambda msg: None)

    def flush_errors(filename: str | None, new_messages: list[str], is_serious: bool) -> None:
        # We never flush errors in the worker, we send them back to coordinator.
        pass

    try:
        return BuildManager(
            data_dir,
            search_paths,
            ignore_prefix=os.getcwd(),
            source_set=source_set,
            reports=None,
            options=options,
            version_id=__version__,
            plugin=plugin,
            plugins_snapshot=snapshot,
            errors=ctx.errors,
            error_formatter=None if options.output is None else OUTPUT_CHOICES.get(options.output),
            flush_errors=flush_errors,
            fscache=ctx.fscache,
            stdout=sys.stdout,
            stderr=sys.stderr,
            parallel_worker=True,
        )
    except CompileError:
        return None


def console_entry() -> None:
    main(sys.argv[1:])
