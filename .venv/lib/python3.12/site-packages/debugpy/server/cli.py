# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License. See LICENSE in the project root
# for license information.

import json
import os
import re
import sys
from importlib.util import find_spec
from typing import Any, Union, Tuple, Dict, Literal

# debugpy.__main__ should have preloaded pydevd properly before importing this module.
# Otherwise, some stdlib modules above might have had imported threading before pydevd
# could perform the necessary detours in it.
assert "pydevd" in sys.modules
import pydevd

# Note: use the one bundled from pydevd so that it's invisible for the user.
from _pydevd_bundle import pydevd_runpy as runpy

import debugpy
import debugpy.server
from debugpy.common import log, sockets
from debugpy.server import api

TargetKind = Literal["file", "module", "code", "pid"]

TARGET = "<filename> | -m <module> | -c <code> | --pid <pid>"

HELP = """debugpy {0}
See https://aka.ms/debugpy for documentation.

Usage: debugpy --listen | --connect
               [<host>:]<port>
               [--wait-for-client]
               [--configure-<name> <value>]...
               [--log-to <path>] [--log-to-stderr]
               [--parent-session-pid <pid>]]
               [--adapter-access-token <token>]
               [--disable-sys-remote-exec]]
               {1}
               [<arg>]...
""".format(
    debugpy.__version__, TARGET
)


# Changes here should be aligned with the public API CliOptions.
class Options(object):
    mode: Union[Literal["connect", "listen"], None] = None
    address: Union[Tuple[str, int], None] = None
    log_to = None
    log_to_stderr = False
    target: Union[str, None] = None
    target_kind: Union[TargetKind, None] = None
    wait_for_client = False
    adapter_access_token = None
    config: Dict[str, Any] = {}
    parent_session_pid: Union[int, None] = None
    disable_sys_remote_exec = False


options = Options()
options.config = {"qt": "none", "subProcess": True}


def in_range(parser, start, stop):
    def parse(s):
        n = parser(s)
        if start is not None and n < start:
            raise ValueError("must be >= {0}".format(start))
        if stop is not None and n >= stop:
            raise ValueError("must be < {0}".format(stop))
        return n

    return parse


pid = in_range(int, 0, None)


def print_help_and_exit(switch, it):
    print(HELP, file=sys.stderr)
    sys.exit(0)


def print_version_and_exit(switch, it):
    print(debugpy.__version__)
    sys.exit(0)


def set_arg(varname, parser=(lambda x: x)):
    def do(arg, it):
        value = parser(next(it))
        setattr(options, varname, value)

    return do


def set_const(varname, value):
    def do(arg, it):
        setattr(options, varname, value)

    return do


def set_address(mode):
    def do(arg, it):
        if options.address is not None:
            raise ValueError("--listen and --connect are mutually exclusive")

        # It's either host:port, or just port.
        value = next(it)
        host, sep, port = value.rpartition(":")
        host = host.strip("[]")
        if not sep:
            host = sockets.get_default_localhost()
            port = value
        try:
            port = int(port)
        except Exception:
            port = -1
        if not (0 <= port < 2**16):
            raise ValueError("invalid port number")

        options.mode = mode
        options.address = (host, port)

    return do


def set_config(arg, it):
    prefix = "--configure-"
    assert arg.startswith(prefix)
    name = arg[len(prefix) :]
    value = next(it)

    if name not in options.config:
        raise ValueError("unknown property {0!r}".format(name))

    expected_type = type(options.config[name])
    try:
        if expected_type is bool:
            value = {"true": True, "false": False}[value.lower()]
        else:
            value = expected_type(value)
    except Exception:
        raise ValueError("{0!r} must be a {1}".format(name, expected_type.__name__))

    options.config[name] = value


def set_target(kind: TargetKind, parser=(lambda x: x), positional=False):
    def do(arg, it):
        options.target_kind = kind
        target = parser(arg if positional else next(it))

        if isinstance(target, bytes):
            # target may be the code, so, try some additional encodings...
            try:
                target = target.decode(sys.getfilesystemencoding())
            except UnicodeDecodeError:
                try:
                    target = target.decode("utf-8")
                except UnicodeDecodeError:
                    import locale

                    target = target.decode(locale.getpreferredencoding(False))

        options.target = target

    return do


# fmt: off
switches = [
    # Switch                    Placeholder         Action
    # ======                    ===========         ======

    # Switches that are documented for use by end users.
    ("-(\\?|h|-help)",          None,               print_help_and_exit),
    ("-(V|-version)",           None,               print_version_and_exit),
    ("--log-to" ,               "<path>",           set_arg("log_to")),
    ("--log-to-stderr",         None,               set_const("log_to_stderr", True)),
    ("--listen",                "<address>",        set_address("listen")),
    ("--connect",               "<address>",        set_address("connect")),
    ("--wait-for-client",       None,               set_const("wait_for_client", True)),
    ("--configure-.+",          "<value>",          set_config),
    ("--parent-session-pid",    "<pid>",            set_arg("parent_session_pid", lambda x: int(x) if x else None)),
    ("--adapter-access-token",   "<token>",         set_arg("adapter_access_token")),
    ("--disable-sys-remote-exec", None,             set_const("disable_sys_remote_exec", True)),

    # Targets. The "" entry corresponds to positional command line arguments,
    # i.e. the ones not preceded by any switch name.
    ("",                        "<filename>",       set_target("file", positional=True)),
    ("-m",                      "<module>",         set_target("module")),
    ("-c",                      "<code>",           set_target("code")),
    ("--pid",                   "<pid>",            set_target("pid", pid)),
]
# fmt: on


# Consume all the args from argv
def consume_argv():
    while len(sys.argv) >= 2:
        value = sys.argv[1]
        del sys.argv[1]
        yield value


# Consume all the args from a given list
def consume_args(args: list):
    if args is sys.argv:
        yield from consume_argv()
    else:
        while args:
            value = args[0]
            del args[0]
            yield value


# Parse the args from the command line, then from the environment.
# Args from the environment are only used if they are not already set from the command line.
def parse_args():

    # keep track of the switches we've seen so far
    seen = set()

    parse_args_from_command_line(seen)
    parse_args_from_environment(seen)

    # if the target is not set, or is empty, this is an error
    if options.target is None or options.target == "":
        raise ValueError("missing target: " + TARGET)

    if options.mode is None:
        raise ValueError("either --listen or --connect is required")
    if options.adapter_access_token is not None and options.mode != "connect":
        raise ValueError("--adapter-access-token requires --connect")
    if options.parent_session_pid is not None and options.mode != "connect":
        raise ValueError("--parent-session-pid requires --connect")
    if options.target_kind == "pid" and options.wait_for_client:
        raise ValueError("--pid does not support --wait-for-client")

    assert options.target_kind is not None
    assert options.address is not None


def parse_args_from_command_line(seen: set):
    parse_args_helper(sys.argv, seen)


def parse_args_from_environment(seenFromCommandLine: set):
    args = os.environ.get("DEBUGPY_EXTRA_ARGV")
    if not args:
        return

    argsList = args.split()

    seenFromEnvironment = set()
    parse_args_helper(argsList, seenFromCommandLine, seenFromEnvironment, True)


def parse_args_helper(
    args: list,
    seenFromCommandLine: set,
    seenFromEnvironment: set = set(),
    isFromEnvironment=False,
):
    iterator = consume_args(args)

    while True:
        try:
            arg = next(iterator)
        except StopIteration:
            break

        switch = arg
        if not switch.startswith("-"):
            switch = ""
        for pattern, placeholder, action in switches:
            if re.match("^(" + pattern + ")$", switch):
                break
        else:
            raise ValueError("unrecognized switch " + switch)

        # if we're parsing from the command line, and we've already seen the switch on the command line, this is an error
        if not isFromEnvironment and switch in seenFromCommandLine:
            raise ValueError("duplicate switch on command line: " + switch)
        # if we're parsing from the environment, and we've already seen the switch in the environment, this is an error
        elif isFromEnvironment and switch in seenFromEnvironment:
            raise ValueError("duplicate switch from environment: " + switch)
        # if we're parsing from the environment, and we've already seen the switch on the command line, skip it, since command line takes precedence
        elif isFromEnvironment and switch in seenFromCommandLine:
            continue
        # otherwise, the switch is new, so add it to the appropriate set
        else:
            if isFromEnvironment:
                seenFromEnvironment.add(switch)
            else:
                seenFromCommandLine.add(switch)

        # process the switch, running the corresponding action
        try:
            action(arg, iterator)
        except StopIteration:
            assert placeholder is not None
            raise ValueError("{0}: missing {1}".format(switch, placeholder))
        except Exception as exc:
            raise ValueError("invalid {0} {1}: {2}".format(switch, placeholder, exc))

        # If we're parsing the command line, we're done after we've processed the target
        # Otherwise, we need to keep parsing until all args are consumed, since the target may be set from the command line
        # already, but there might be additional args in the environment that we want to process.
        if not isFromEnvironment and options.target is not None:
            break


def start_debugging(argv_0):
    # We need to set up sys.argv[0] before invoking either listen() or connect(),
    # because they use it to report the "process" event. Thus, we can't rely on
    # run_path() and run_module() doing that, even though they will eventually.
    sys.argv[0] = argv_0

    log.debug("sys.argv after patching: {0!r}", sys.argv)

    debugpy.configure(options.config)

    if os.environ.get("DEBUGPY_RUNNING", "false") != "true":
        if options.mode == "listen" and options.address is not None:
            debugpy.listen(options.address)
        elif options.mode == "connect" and options.address is not None:
            debugpy.connect(options.address, access_token=options.adapter_access_token, parent_session_pid=options.parent_session_pid)
        else:
            raise AssertionError(repr(options.mode))

        if options.wait_for_client:
            debugpy.wait_for_client()

    os.environ["DEBUGPY_RUNNING"] = "true"


def run_file():
    target = options.target
    start_debugging(target)

    # run_path has one difference with invoking Python from command-line:
    # if the target is a file (rather than a directory), it does not add its
    # parent directory to sys.path. Thus, importing other modules from the
    # same directory is broken unless sys.path is patched here.

    if target is not None and os.path.isfile(target):
        dir = os.path.dirname(target)
        sys.path.insert(0, dir)
    else:
        log.debug("Not a file: {0!r}", target)

    log.describe_environment("Pre-launch environment:")

    log.info("Running file {0!r}", target)
    runpy.run_path(target, run_name="__main__")


def run_module():
    # Add current directory to path, like Python itself does for -m. This must
    # be in place before trying to use find_spec below to resolve submodules.
    sys.path.insert(0, str(""))

    # We want to do the same thing that run_module() would do here, without
    # actually invoking it.
    argv_0 = sys.argv[0]
    try:
        spec = None if options.target is None else find_spec(options.target)
        if spec is not None:
            argv_0 = spec.origin
    except Exception:
        log.swallow_exception("Error determining module path for sys.argv")

    start_debugging(argv_0)
    log.describe_environment("Pre-launch environment:")
    log.info("Running module {0!r}", options.target)

    # Docs say that runpy.run_module is equivalent to -m, but it's not actually
    # the case for packages - -m sets __name__ to "__main__", but run_module sets
    # it to "pkg.__main__". This breaks everything that uses the standard pattern
    # __name__ == "__main__" to detect being run as a CLI app. On the other hand,
    # runpy._run_module_as_main is a private function that actually implements -m.
    try:
        run_module_as_main = runpy._run_module_as_main
    except AttributeError:
        log.warning("runpy._run_module_as_main is missing, falling back to run_module.")
        runpy.run_module(options.target, alter_sys=True)
    else:
        run_module_as_main(options.target, alter_argv=True)


def run_code():
    if options.target is not None:
        # Add current directory to path, like Python itself does for -c.
        sys.path.insert(0, str(""))
        code = compile(options.target, str("<string>"), str("exec"))

        start_debugging(str("-c"))

        log.describe_environment("Pre-launch environment:")
        log.info("Running code:\n\n{0}", options.target)

        eval(code, {})
    else:
        log.error("No target to run.")


def attach_to_pid():
    pid = options.target
    log.info("Attaching to process with PID={0}", pid)

    encode = lambda s: list(bytearray(s.encode("utf-8"))) if s is not None else None

    script_dir = os.path.dirname(debugpy.server.__file__)
    assert os.path.exists(script_dir)
    script_dir = encode(script_dir)

    setup = {
        "mode": options.mode,
        "address": options.address,
        "wait_for_client": options.wait_for_client,
        "log_to": options.log_to,
        "adapter_access_token": options.adapter_access_token,
    }
    setup = encode(json.dumps(setup))

    python_code = """
import codecs;
import json;
import sys;

decode = lambda s: codecs.utf_8_decode(bytearray(s))[0] if s is not None else None;

script_dir = decode({script_dir});
setup = json.loads(decode({setup}));

sys.path.insert(0, script_dir);
import attach_pid_injected;
del sys.path[0];

attach_pid_injected.attach(setup);
"""
    python_code = (
        python_code.replace("\r", "")
        .replace("\n", "")
        .format(script_dir=script_dir, setup=setup)
    )

    # attempt pep 768 style code injection
    if (not options.disable_sys_remote_exec) and hasattr(sys, "remote_exec"):
        tmp_file_path = ""
        try:
            import tempfile

            with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
                tmp_file_path = tmp_file.name
                log.info(
                    "Attempting to inject code at '{tmp_file_path}' using sys.remote_exec()",
                    tmp_file_path=tmp_file_path,
                )
                tmp_file.write(python_code.encode())
                tmp_file.write(
                    """import os;os.remove("{tmp_file_path}");""".format(
                        tmp_file_path=tmp_file_path
                    ).encode()
                )
                tmp_file.flush()
                tmp_file.close()
                sys.remote_exec(pid, tmp_file_path)
            return
        except Exception as e:
            if os.path.exists(tmp_file_path):
                os.remove(tmp_file_path)
            log.warning(
                'Injecting code using sys.remote_exec() failed with error:\n"{e}"\nWill reattempt using pydevd.\n',
                e=e,
            )

    log.info("Code to be injected: \n{0}", python_code.replace(";", ";\n"))

    # pydevd restriction on characters in injected code.
    assert not (
        {'"', "'", "\r", "\n"} & set(python_code)
    ), "Injected code should not contain any single quotes, double quotes, or newlines."

    pydevd_attach_to_process_path = os.path.join(
        os.path.dirname(pydevd.__file__), "pydevd_attach_to_process"
    )

    assert os.path.exists(pydevd_attach_to_process_path)
    sys.path.append(pydevd_attach_to_process_path)

    try:
        import add_code_to_python_process  # noqa

        log.info("Injecting code into process with PID={0} ...", pid)
        add_code_to_python_process.run_python_code(
            pid,
            python_code,
            connect_debugger_tracing=True,
            show_debug_info=int(os.getenv("DEBUGPY_ATTACH_BY_PID_DEBUG_INFO", "0")),
        )
    except Exception:
        log.reraise_exception("Code injection into PID={0} failed:", pid)
    log.info("Code injection into PID={0} completed.", pid)


def main():
    original_argv = list(sys.argv)
    try:
        parse_args()
    except Exception as exc:
        print(str(HELP) + str("\nError: ") + str(exc), file=sys.stderr)
        sys.exit(2)

    if options.log_to is not None:
        debugpy.log_to(options.log_to)
    if options.log_to_stderr:
        debugpy.log_to(sys.stderr)

    api.ensure_logging()

    log.info(
        str("sys.argv before parsing: {0!r}\n" "         after parsing:  {1!r}"),
        original_argv,
        sys.argv,
    )

    try:
        if options.target_kind is not None:
            run = {
                "file": run_file,
                "module": run_module,
                "code": run_code,
                "pid": attach_to_pid,
            }[options.target_kind]
            run()
    except SystemExit as exc:
        log.reraise_exception(
            "Debuggee exited via SystemExit: {0!r}", exc.code, level="debug"
        )
