import os
import sys

from llvmlite import ir

from numba.core import types, utils, config, cgutils, errors
from numba import gdb, gdb_init, gdb_breakpoint
from numba.core.extending import overload, intrinsic

_path = os.path.dirname(__file__)

_platform = sys.platform
_unix_like = (_platform.startswith('linux') or
              _platform.startswith('darwin') or
              ('bsd' in _platform))


def _confirm_gdb(need_ptrace_attach=True):
    """
    Set need_ptrace_attach to True/False to indicate whether the ptrace attach
    permission is needed for this gdb use case. Mode 0 (classic) or 1
    (restricted ptrace) is required if need_ptrace_attach is True. See:
    https://www.kernel.org/doc/Documentation/admin-guide/LSM/Yama.rst
    for details on the modes.
    """
    if not _unix_like:
        msg = 'gdb support is only available on unix-like systems'
        raise errors.NumbaRuntimeError(msg)
    gdbloc = config.GDB_BINARY
    if not (os.path.exists(gdbloc) and os.path.isfile(gdbloc)):
        msg = ('Is gdb present? Location specified (%s) does not exist. The gdb'
               ' binary location can be set using Numba configuration, see: '
               'https://numba.readthedocs.io/en/stable/reference/envvars.html'  # noqa: E501
               )
        raise RuntimeError(msg % config.GDB_BINARY)
    # Is Yama being used as a kernel security module and if so is ptrace_scope
    # limited? In this case ptracing non-child processes requires special
    # permission so raise an exception.
    ptrace_scope_file = os.path.join(os.sep, 'proc', 'sys', 'kernel', 'yama',
                                     'ptrace_scope')
    has_ptrace_scope = os.path.exists(ptrace_scope_file)
    if has_ptrace_scope:
        with open(ptrace_scope_file, 'rt') as f:
            value = f.readline().strip()
        if need_ptrace_attach and value not in ("0", "1"):
            msg = ("gdb can launch but cannot attach to the executing program"
                   " because ptrace permissions have been restricted at the "
                   "system level by the Linux security module 'Yama'.\n\n"
                   "Documentation for this module and the security "
                   "implications of making changes to its behaviour can be "
                   "found in the Linux Kernel documentation "
                   "https://www.kernel.org/doc/Documentation/admin-guide/LSM/Yama.rst"    # noqa: E501
                   "\n\nDocumentation on how to adjust the behaviour of Yama "
                   "on Ubuntu Linux with regards to 'ptrace_scope' can be "
                   "found here "
                   "https://wiki.ubuntu.com/Security/Features#ptrace.")
            raise RuntimeError(msg)


@overload(gdb)
def hook_gdb(*args):
    _confirm_gdb()
    gdbimpl = gen_gdb_impl(args, True)

    def impl(*args):
        gdbimpl()
    return impl


@overload(gdb_init)
def hook_gdb_init(*args):
    _confirm_gdb()
    gdbimpl = gen_gdb_impl(args, False)

    def impl(*args):
        gdbimpl()
    return impl


def init_gdb_codegen(cgctx, builder, signature, args,
                     const_args, do_break=False):

    int8_t = ir.IntType(8)
    int32_t = ir.IntType(32)
    intp_t = ir.IntType(utils.MACHINE_BITS)
    char_ptr = ir.PointerType(ir.IntType(8))
    zero_i32t = int32_t(0)

    mod = builder.module
    pid = cgutils.alloca_once(builder, int32_t, size=1)

    # 32bit pid, 11 char max + terminator
    pidstr = cgutils.alloca_once(builder, int8_t, size=12)

    # str consts
    intfmt = cgctx.insert_const_string(mod, '%d')
    gdb_str = cgctx.insert_const_string(mod, config.GDB_BINARY)
    attach_str = cgctx.insert_const_string(mod, 'attach')

    new_args = []
    # add break point command to known location
    # this command file thing is due to commands attached to a breakpoint
    # requiring an interactive prompt
    # https://sourceware.org/bugzilla/show_bug.cgi?id=10079
    new_args.extend(['-x', os.path.join(_path, 'cmdlang.gdb')])
    # issue command to continue execution from sleep function
    new_args.extend(['-ex', 'c'])
    # then run the user defined args if any
    if any([not isinstance(x, types.StringLiteral) for x in const_args]):
        raise errors.RequireLiteralValue(const_args)
    new_args.extend([x.literal_value for x in const_args])
    cmdlang = [cgctx.insert_const_string(mod, x) for x in new_args]

    # insert getpid, getpid is always successful, call without concern!
    fnty = ir.FunctionType(int32_t, tuple())
    getpid = cgutils.get_or_insert_function(mod, fnty, "getpid")

    # insert snprintf
    # int snprintf(char *str, size_t size, const char *format, ...);
    fnty = ir.FunctionType(
        int32_t, (char_ptr, intp_t, char_ptr), var_arg=True)
    snprintf = cgutils.get_or_insert_function(mod, fnty, "snprintf")

    # insert fork
    fnty = ir.FunctionType(int32_t, tuple())
    fork = cgutils.get_or_insert_function(mod, fnty, "fork")

    # insert execl
    fnty = ir.FunctionType(int32_t, (char_ptr, char_ptr), var_arg=True)
    execl = cgutils.get_or_insert_function(mod, fnty, "execl")

    # insert sleep
    fnty = ir.FunctionType(int32_t, (int32_t,))
    sleep = cgutils.get_or_insert_function(mod, fnty, "sleep")

    # insert break point
    fnty = ir.FunctionType(ir.VoidType(), tuple())
    breakpoint = cgutils.get_or_insert_function(mod, fnty,
                                                "numba_gdb_breakpoint")

    # do the work
    parent_pid = builder.call(getpid, tuple())
    builder.store(parent_pid, pid)
    pidstr_ptr = builder.gep(pidstr, [zero_i32t], inbounds=True)
    pid_val = builder.load(pid)

    # call snprintf to write the pid into a char *
    stat = builder.call(
        snprintf, (pidstr_ptr, intp_t(12), intfmt, pid_val))
    invalid_write = builder.icmp_signed('>', stat, int32_t(12))
    with builder.if_then(invalid_write, likely=False):
        msg = "Internal error: `snprintf` buffer would have overflowed."
        cgctx.call_conv.return_user_exc(builder, RuntimeError, (msg,))

    # fork, check pids etc
    child_pid = builder.call(fork, tuple())
    fork_failed = builder.icmp_signed('==', child_pid, int32_t(-1))
    with builder.if_then(fork_failed, likely=False):
        msg = "Internal error: `fork` failed."
        cgctx.call_conv.return_user_exc(builder, RuntimeError, (msg,))

    is_child = builder.icmp_signed('==', child_pid, zero_i32t)
    with builder.if_else(is_child) as (then, orelse):
        with then:
            # is child
            nullptr = ir.Constant(char_ptr, None)
            gdb_str_ptr = builder.gep(
                gdb_str, [zero_i32t], inbounds=True)
            attach_str_ptr = builder.gep(
                attach_str, [zero_i32t], inbounds=True)
            cgutils.printf(
                builder, "Attaching to PID: %s\n", pidstr)
            buf = (
                gdb_str_ptr,
                gdb_str_ptr,
                attach_str_ptr,
                pidstr_ptr)
            buf = buf + tuple(cmdlang) + (nullptr,)
            builder.call(execl, buf)
        with orelse:
            # is parent
            builder.call(sleep, (int32_t(10),))
            # if breaking is desired, break now
            if do_break is True:
                builder.call(breakpoint, tuple())


def gen_gdb_impl(const_args, do_break):
    @intrinsic
    def gdb_internal(tyctx):
        function_sig = types.void()

        def codegen(cgctx, builder, signature, args):
            init_gdb_codegen(cgctx, builder, signature, args, const_args,
                             do_break=do_break)
            return cgctx.get_constant(types.none, None)
        return function_sig, codegen
    return gdb_internal


@overload(gdb_breakpoint)
def hook_gdb_breakpoint():
    """
    Adds the Numba break point into the source
    """
    if not sys.platform.startswith('linux'):
        raise RuntimeError('gdb is only available on linux')
    bp_impl = gen_bp_impl()

    def impl():
        bp_impl()
    return impl


def gen_bp_impl():
    @intrinsic
    def bp_internal(tyctx):
        function_sig = types.void()

        def codegen(cgctx, builder, signature, args):
            mod = builder.module
            fnty = ir.FunctionType(ir.VoidType(), tuple())
            breakpoint = cgutils.get_or_insert_function(mod, fnty,
                                                        "numba_gdb_breakpoint")
            builder.call(breakpoint, tuple())
            return cgctx.get_constant(types.none, None)
        return function_sig, codegen
    return bp_internal
