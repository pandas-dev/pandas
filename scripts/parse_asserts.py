#!/usr/bin/env python

import re
import os
import fnmatch
import ast
import argparse
import inspect
import sys
import tempfile
import subprocess
import operator

try:
    from importlib import import_module
except ImportError:
    import_module = __import__


from numpy import nan as NA
from pandas import DataFrame
from pandas.core.config import option_context


def parse_interp_string(node):
    assert isinstance(node, ast.BinOp)
    assert isinstance(node.op, ast.Mod)
    assert isinstance(node.left, ast.Str)
    return node.left.s


def parse_format_string(node):
    assert isinstance(node, ast.Call)
    assert isinstance(node.func, ast.Attribute)
    assert isinstance(node.func.value, ast.Str)
    return node.func.value.s


def try_parse_raise_arg(node):
    try:
        # string
        v = node.s
    except AttributeError:
        try:
            # interpolated string
            v = parse_interp_string(node)
        except AssertionError:
            try:
                # format spec string
                v = parse_format_string(node)
            except AssertionError:
                # otherwise forget it (general expr node)
                v = node
    return v


def parse_file(pyfile, asserts):
    with open(pyfile, 'r') as pyf:
        source = pyf.read()

    try:
        parsed = ast.parse(source, pyfile, 'exec')
    except SyntaxError:
        return

    for node in ast.walk(parsed):
        if isinstance(node, ast.Raise):
            k = pyfile, node.lineno, node.col_offset

            try:
                # try to get the name of the exception constructor
                asserts[k] = [node.type.func.id]
            except AttributeError:
                # not a constructor
                asserts[k] = [NA]
            else:
                # is constructor, try parsing its contents
                try:
                    # function arguments
                    args = node.type.args

                    try:
                        # try to get the first argument
                        arg = args[0]
                        v = try_parse_raise_arg(arg)
                        asserts[k].append(v)
                    except IndexError:
                        # no arguments (e.g., raise Exception())
                        asserts[k].append(NA)

                except AttributeError:
                    # no arguments (e.g., raise Exception)
                    asserts[k].append(NA)


def path_matches(path, pattern):
    return re.search(pattern, path) is not None


def regex_or(*patterns):
    return '({0})'.format('|'.join(patterns))


def get_asserts_from_path(path, file_filters, dir_filters):
    if file_filters is None:
        file_filters = 'test', '__init__.py'

    file_filters = regex_or(*file_filters)

    if dir_filters is None:
        dir_filters = 'build', '.tox', 'test', '.*\.egg.*'

    dir_filters = regex_or(*dir_filters)

    asserts = {}

    if os.path.isfile(path):
        parse_file(path, asserts)
        return asserts

    for root, _, filenames in os.walk(path):
        full_names = []

        if not path_matches(root, dir_filters):
            full_names = [os.path.join(root, fn) for fn in filenames
                          if not path_matches(fn, file_filters)]

        if full_names:
            pyfiles = fnmatch.filter(full_names, '*.py')

            if pyfiles:
                for pyfile in pyfiles:
                    parse_file(pyfile, asserts)

    return asserts


def obj_path_from_string(dotted_name, full_path):
    try:
        obj = import_module(dotted_name)
    except ImportError:
        splits_ville = dotted_name.split('.')
        module_name, obj_name = splits_ville[:-1], splits_ville[-1]
        module_name = '.'.join(module_name)

        try:
            module = import_module(module_name)
        except ImportError:
            raise ImportError("'{0}' is not a valid Python "
                              "module".format(module_name))
        else:
            try:
                obj = getattr(module, obj_name)
            except AttributeError:
                raise AttributeError("")

    if full_path:
        path = inspect.getabsfile(obj)
    else:
        path = inspect.getfile(obj)

    if path.endswith('pyc'):
        path = path.strip('c')
    return os.path.dirname(path)


def get_asserts_from_obj(dotted_name, file_filters, dir_filters, full_path):
    path = obj_path_from_string(dotted_name, full_path)
    return get_asserts_from_path(path, file_filters, dir_filters)


def asserts_to_frame(asserts):
    index, values = zip(*asserts.iteritems())
    values = map(lambda x: list(reduce(operator.concat, map(list, x))),
                 asserts.iteritems())
    columns = 'filename', 'line', 'col', 'code', 'msg'
    df = DataFrame(values, columns=columns).fillna(NA)
    return df


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--type', default='all',
                        choices=('all', 'a', 'empty', 'e', 'nonempty', 'n'),
                        help='The type of nodes you want to look for')
    parser.add_argument('-m', '--module', default='pandas',
                        help=('The name of a module or file to search for '
                              'nodes in'))
    parser.add_argument('-i', '--file-filters', default=None, nargs='*',
                        help=("A list of regular expressions describing files "
                              "you want to ignore"))
    parser.add_argument('-d', '--dir-filters', default=None, nargs='*',
                        help=('A list of regular expressions describing'
                              ' directories you want to ignore'))
    parser.add_argument('-s', '--sparse-filename', action='store_true',
                        help=('Use multi_sparse = False to show the '
                              'resulting DataFrame'))
    parser.add_argument('-p', '--full-path', action='store_true',
                        help=('Display the entire path of the file if this '
                              'is given'))
    parser.add_argument('-k', '--exception-types', nargs='*',
                        help='The types of exceptions to report')
    parser.add_argument('-b', '--sort-by', default='line', nargs='*',
                        help=('A list of columns or index levels you want to '
                              'sort by'))
    return parser.parse_args()


def _build_exc_regex(exc_list):
    return r'(.*(?:{0}).*)'.format('|'.join(exc_list))


def main(args):
    asserts = get_asserts_from_obj(args.module, args.file_filters,
                                   args.dir_filters, args.full_path)

    if not asserts:
        print "No asserts found in '{0}'".format(args.module)
        return 0

    df = asserts_to_frame(asserts)

    try:
        df.sortlevel(args.sort_by, inplace=True)
    except Exception:
        df.sort(args.sort_by, inplace=True)

    atype = args.type

    msg = 'No'

    if atype.startswith('e'):
        ind = df.msg.isnull()
        msg += ' empty'
    elif atype.startswith('n'):
        ind = df.msg.notnull()
        msg += ' nonempty'
    else:
        ind = slice(None)

    df = df[ind]
    df.sort_index(inplace=True)

    if df.empty:
        print "{0} {1} found in '{2}'".format(msg, args.exception_types,
                                              args.module)
        return 0
    max_cols = int(df.msg.map(lambda x: len(repr(x))).max())
    with option_context('display.multi_sparse', args.sparse_filename,
                        'display.max_colwidth', max_cols,
                        'display.max_seq_items', max_cols):
        if args.exception_types is not None:
            regex = _build_exc_regex(args.exception_types)
            vals = df.code.str.match(regex, re.I)
            df = df[vals.str[0].notnull()]

            if df.empty:
                msg = "{0} {1} found in '{2}'".format(msg,
                                                      args.exception_types,
                                                      args.module)
                print msg
                return 0

            with tempfile.NamedTemporaryFile() as tmpf:
                df.to_string(buf=tmpf)
                return subprocess.call([os.environ.get('PAGER', 'less'),
                                        tmpf.name])
    return df


if __name__ == '__main__':
    sys.exit(main(parse_args()))
