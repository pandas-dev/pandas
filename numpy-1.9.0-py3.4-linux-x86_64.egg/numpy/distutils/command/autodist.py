"""This module implements additional tests ala autoconf which can be useful.

"""
from __future__ import division, absolute_import, print_function


# We put them here since they could be easily reused outside numpy.distutils

def check_inline(cmd):
    """Return the inline identifier (may be empty)."""
    cmd._check_compiler()
    body = """
#ifndef __cplusplus
static %(inline)s int static_func (void)
{
    return 0;
}
%(inline)s int nostatic_func (void)
{
    return 0;
}
#endif"""

    for kw in ['inline', '__inline__', '__inline']:
        st = cmd.try_compile(body % {'inline': kw}, None, None)
        if st:
            return kw

    return ''

def check_compiler_gcc4(cmd):
    """Return True if the C compiler is GCC 4.x."""
    cmd._check_compiler()
    body = """
int
main()
{
#if (! defined __GNUC__) || (__GNUC__ < 4)
#error gcc >= 4 required
#endif
}
"""
    return cmd.try_compile(body, None, None)
