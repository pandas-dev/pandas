""" Based on numpy's approach to exposing compiler features via a config header.
Unfortunately that file is not exposed, so re-implement the portions we need.
"""
import os
import textwrap

OPTIONAL_FUNCTION_ATTRIBUTES = [
    ("HAVE_ATTRIBUTE_OPTIMIZE_OPT_3", '__attribute__((optimize("O3")))')
]


def _get_compiler_list(cmd):
    """ Return the compiler command as a list of strings. Distutils provides a
    wildly inconsistent API here:
      - UnixCCompiler returns a list
      - MSVCCompiler intentionally doesn't set this variable
      - CygwinCompiler returns a string

    As we are focused on identifying gcc vs clang right now, we ignore MSVC's
    bad result and convert all results into lists of strings
    """
    compiler = getattr(cmd.compiler, "compiler", "")
    if isinstance(compiler, str):
        compiler = compiler.split()
    return compiler


def is_gcc(cmd):
    return any("gcc" in x for x in _get_compiler_list(cmd))


def is_clang(cmd):
    return any("clang" in x for x in _get_compiler_list(cmd))


def check_inline(cmd):
    """Return the inline identifier (may be empty)."""
    cmd._check_compiler()
    body = textwrap.dedent(
        """
        #ifndef __cplusplus
        static %(inline)s int static_func (void)
        {
            return 0;
        }
        %(inline)s int nostatic_func (void)
        {
            return 0;
        }
        #endif
        int main(void) {
            int r1 = static_func();
            int r2 = nostatic_func();
            return r1 + r2;
        }
        """
    )

    for kw in ["inline", "__inline__", "__inline"]:
        st = cmd.try_compile(body % {"inline": kw}, None, None)
        if st:
            return kw

    return ""


def check_gcc_function_attribute(cmd, attribute, name):
    """Return True if the given function attribute is supported."""
    cmd._check_compiler()
    if is_gcc(cmd):
        pragma = '#pragma GCC diagnostic error "-Wattributes"'
    elif is_clang(cmd):
        pragma = '#pragma clang diagnostic error "-Wattributes"'
    else:
        pragma = ""

    body = (
        textwrap.dedent(
            """
        %s

        int %s %s(void*);

        int main(void)
        {
            return 0;
        }
        """
        )
        % (pragma, attribute, name)
    )
    return cmd.try_compile(body, None, None) != 0


def create_config_h(config):
    dirname = os.path.dirname(__file__)
    config_h = os.path.join(dirname, "bn_config.h")

    if (
        os.path.exists(config_h)
        and os.stat(__file__).st_mtime < os.stat(config_h).st_mtime
    ):
        return

    output = []

    for config_attr, func_attr in OPTIONAL_FUNCTION_ATTRIBUTES:
        if check_gcc_function_attribute(config, func_attr, config_attr.lower()):
            output.append((config_attr, "1"))
        else:
            output.append((config_attr, "0"))

    inline_alias = check_inline(config)

    with open(config_h, "w") as f:
        for setting in output:
            f.write("#define {} {}\n".format(*setting))

        if inline_alias == "inline":
            f.write("/* undef inline */\n")
        else:
            f.write("#define inline {}\n".format(inline_alias))
