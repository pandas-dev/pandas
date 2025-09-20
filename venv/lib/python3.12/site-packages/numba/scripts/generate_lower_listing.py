"""
Generate documentation for all registered implementation for lowering
using reStructured text.
"""


from subprocess import check_output

import os.path
try:
    from StringIO import StringIO       # py2
except ImportError:
    from io import StringIO
from collections import defaultdict
import inspect
from functools import partial

import numba
from numba.core.registry import cpu_target


def git_hash():
    out = check_output(['git', 'log', "--pretty=format:'%H'", '-n', '1'])
    return out.decode('ascii').strip("'\"")


def get_func_name(fn):
    return getattr(fn, '__qualname__', fn.__name__)


def gather_function_info(backend):
    fninfos = defaultdict(list)
    basepath = os.path.dirname(os.path.dirname(numba.__file__))
    for fn, osel in backend._defns.items():
        for sig, impl in osel.versions:
            info = {}
            fninfos[fn].append(info)
            info['fn'] = fn
            info['sig'] = sig
            code, firstlineno = inspect.getsourcelines(impl)
            path = inspect.getsourcefile(impl)
            info['impl'] = {
                'name': get_func_name(impl),
                'filename': os.path.relpath(path, start=basepath),
                'lines': (firstlineno, firstlineno + len(code) - 1),
                'docstring': impl.__doc__
            }

    return fninfos


def bind_file_to_print(fobj):
    return partial(print, file=fobj)


def format_signature(sig):
    def fmt(c):
        try:
            return c.__name__
        except AttributeError:
            return repr(c).strip('\'"')
    out = tuple(map(fmt, sig))
    return '`({0})`'.format(', '.join(out))


github_url = ('https://github.com/numba/numba/blob/'
              '{commit}/{path}#L{firstline}-L{lastline}')

description = """
This lists all lowering definition registered to the CPU target.
Each subsection corresponds to a Python function that is supported by numba
nopython mode. These functions have one or more lower implementation with
different signatures. The compiler chooses the most specific implementation
from all overloads.
"""


def format_function_infos(fninfos):
    buf = StringIO()
    try:
        print = bind_file_to_print(buf)

        title_line = "Lowering Listing"
        print(title_line)
        print('=' * len(title_line))

        print(description)

        commit = git_hash()

        def format_fname(fn):
            try:
                fname = "{0}.{1}".format(fn.__module__, get_func_name(fn))
            except AttributeError:
                fname = repr(fn)
            return fn, fname

        for fn, fname in sorted(map(format_fname, fninfos), key=lambda x: x[1]):
            impinfos = fninfos[fn]
            header_line = "``{0}``".format(fname)
            print(header_line)
            print('-' * len(header_line))
            print()

            formatted_sigs = map(
                lambda x: format_signature(x['sig']), impinfos)
            sorted_impinfos = sorted(zip(formatted_sigs, impinfos),
                                     key=lambda x: x[0])

            col_signatures = ['Signature']
            col_urls = ['Definition']

            for fmtsig, info in sorted_impinfos:
                impl = info['impl']

                filename = impl['filename']
                lines = impl['lines']
                fname = impl['name']

                source = '{0} lines {1}-{2}'.format(filename, *lines)
                link = github_url.format(commit=commit, path=filename,
                                         firstline=lines[0], lastline=lines[1])
                url = '``{0}`` `{1} <{2}>`_'.format(fname, source, link)

                col_signatures.append(fmtsig)
                col_urls.append(url)

            # table formatting
            max_width_col_sig = max(map(len, col_signatures))
            max_width_col_url = max(map(len, col_urls))
            padding = 2
            width_col_sig = padding * 2 + max_width_col_sig
            width_col_url = padding * 2 + max_width_col_url
            line_format = "{{0:^{0}}}  {{1:^{1}}}".format(width_col_sig,
                                                          width_col_url)
            print(line_format.format('=' * width_col_sig, '=' * width_col_url))
            print(line_format.format(col_signatures[0], col_urls[0]))
            print(line_format.format('=' * width_col_sig, '=' * width_col_url))
            for sig, url in zip(col_signatures[1:], col_urls[1:]):
                print(line_format.format(sig, url))
            print(line_format.format('=' * width_col_sig, '=' * width_col_url))
            print()

        return buf.getvalue()
    finally:
        buf.close()


# Main routine for this module:

def gen_lower_listing(path=None):
    """
    Generate lowering listing to ``path`` or (if None) to stdout.
    """
    cpu_backend = cpu_target.target_context
    cpu_backend.refresh()

    fninfos = gather_function_info(cpu_backend)
    out = format_function_infos(fninfos)

    if path is None:
        print(out)
    else:
        with open(path, 'w') as fobj:
            print(out, file=fobj)


if __name__ == '__main__':
    gen_lower_listing()
