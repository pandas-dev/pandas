#!/usr/bin/env python
"""
Analyze docstrings to detect errors.

If no argument is provided, it does a quick check of docstrings and returns
a csv with all API functions and results of basic checks.

If a function or method is provided in the form "pandas.function",
"pandas.module.class.method", etc. a list of all errors in the docstring for
the specified function or method.

Usage::
    $ ./validate_docstrings.py
    $ ./validate_docstrings.py pandas.DataFrame.head
"""
import os
import sys
import csv
import re
import functools
import collections
import argparse
import pydoc
import inspect
import importlib
import doctest
try:
    from io import StringIO
except ImportError:
    from cStringIO import StringIO
import numpy

BASE_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0, os.path.join(BASE_PATH))
import pandas
from pandas.compat import signature

sys.path.insert(1, os.path.join(BASE_PATH, 'doc', 'sphinxext'))
from numpydoc.docscrape import NumpyDocString
from pandas.io.formats.printing import pprint_thing


PRIVATE_CLASSES = ['NDFrame', 'IndexOpsMixin']


def _load_obj(obj_name):
    for maxsplit in range(1, obj_name.count('.') + 1):
        # TODO when py3 only replace by: module, *func_parts = ...
        func_name_split = obj_name.rsplit('.', maxsplit)
        module = func_name_split[0]
        func_parts = func_name_split[1:]
        try:
            obj = importlib.import_module(module)
        except ImportError:
            pass
        else:
            continue

    if 'module' not in locals():
        raise ImportError('No module can be imported '
                          'from "{}"'.format(obj_name))

    for part in func_parts:
        obj = getattr(obj, part)
    return obj


def _to_original_callable(obj):
    while True:
        if inspect.isfunction(obj) or inspect.isclass(obj):
            f = inspect.getfile(obj)
            if f.startswith('<') and f.endswith('>'):
                return None
            return obj
        if inspect.ismethod(obj):
            obj = obj.__func__
        elif isinstance(obj, functools.partial):
            obj = obj.func
        elif isinstance(obj, property):
            obj = obj.fget
        else:
            return None


def _output_header(title, width=80, char='#'):
    full_line = char * width
    side_len = (width - len(title) - 2) // 2
    adj = '' if len(title) % 2 == 0 else ' '
    title_line = '{side} {title}{adj} {side}'.format(side=char * side_len,
                                                     title=title,
                                                     adj=adj)

    return '\n{full_line}\n{title_line}\n{full_line}\n\n'.format(
        full_line=full_line, title_line=title_line)


class Docstring(object):
    def __init__(self, method_name, method_obj):
        self.method_name = method_name
        self.method_obj = method_obj
        self.raw_doc = method_obj.__doc__ or ''
        self.clean_doc = pydoc.getdoc(self.method_obj)
        self.doc = NumpyDocString(self.clean_doc)

    def __len__(self):
        return len(self.raw_doc)

    @property
    def is_function_or_method(self):
        # TODO(py27): remove ismethod
        return (inspect.isfunction(self.method_obj)
                or inspect.ismethod(self.method_obj))

    @property
    def source_file_name(self):
        fname = inspect.getsourcefile(self.method_obj)
        if fname:
            fname = os.path.relpath(fname, BASE_PATH)
            return fname

    @property
    def source_file_def_line(self):
        try:
            return inspect.getsourcelines(self.method_obj)[-1]
        except OSError:
            pass

    @property
    def github_url(self):
        url = 'https://github.com/pandas-dev/pandas/blob/master/'
        url += '{}#L{}'.format(self.source_file_name,
                               self.source_file_def_line)
        return url

    @property
    def start_blank_lines(self):
        i = None
        if self.raw_doc:
            for i, row in enumerate(self.raw_doc.split('\n')):
                if row.strip():
                    break
        return i

    @property
    def end_blank_lines(self):
        i = None
        if self.raw_doc:
            for i, row in enumerate(reversed(self.raw_doc.split('\n'))):
                if row.strip():
                    break
        return i

    @property
    def double_blank_lines(self):
        prev = True
        for row in self.raw_doc.split('\n'):
            if not prev and not row.strip():
                return True
            prev = row.strip()
        return False

    @property
    def summary(self):
        if not self.doc['Extended Summary'] and len(self.doc['Summary']) > 1:
            return ''
        return ' '.join(self.doc['Summary'])

    @property
    def extended_summary(self):
        if not self.doc['Extended Summary'] and len(self.doc['Summary']) > 1:
            return ' '.join(self.doc['Summary'])
        return ' '.join(self.doc['Extended Summary'])

    @property
    def needs_summary(self):
        return not (bool(self.summary) and bool(self.extended_summary))

    @property
    def doc_parameters(self):
        return collections.OrderedDict((name, (type_, ''.join(desc)))
                                       for name, type_, desc
                                       in self.doc['Parameters'])

    @property
    def signature_parameters(self):
        if inspect.isclass(self.method_obj):
            if (hasattr(self.method_obj, '_accessors')
                and self.method_name.split('.')[-1] in
                self.method_obj._accessors):
                # accessor classes have a signature but don't want to show this
                return tuple()
        try:
            sig = signature(self.method_obj)
        except (TypeError, ValueError):
            # Some objects, mainly in C extensions do not support introspection
            # of the signature
            return tuple()
        params = sig.args
        if sig.varargs:
            params.append("*" + sig.varargs)
        if sig.keywords:
            params.append("**" + sig.keywords)
        params = tuple(params)
        if params and params[0] in ('self', 'cls'):
            return params[1:]
        return params

    @property
    def parameter_mismatches(self):
        errs = []
        signature_params = self.signature_parameters
        doc_params = tuple(self.doc_parameters)
        missing = set(signature_params) - set(doc_params)
        if missing:
            errs.append(
                'Parameters {} not documented'.format(pprint_thing(missing)))
        extra = set(doc_params) - set(signature_params)
        if extra:
            errs.append('Unknown parameters {}'.format(pprint_thing(extra)))
        if (not missing and not extra and signature_params != doc_params
                and not (not signature_params and not doc_params)):
            errs.append('Wrong parameters order. ' +
                        'Actual: {!r}. '.format(signature_params) +
                        'Documented: {!r}'.format(doc_params))

        return errs

    @property
    def correct_parameters(self):
        return not bool(self.parameter_mismatches)

    def parameter_type(self, param):
        return self.doc_parameters[param][0]

    def parameter_desc(self, param):
        return self.doc_parameters[param][1]

    @property
    def see_also(self):
        return collections.OrderedDict((name, ''.join(desc))
                                       for name, desc, _
                                       in self.doc['See Also'])

    @property
    def examples(self):
        return self.doc['Examples']

    @property
    def returns(self):
        return self.doc['Returns']

    @property
    def yields(self):
        return self.doc['Yields']

    @property
    def method_source(self):
        return inspect.getsource(self.method_obj)

    @property
    def first_line_ends_in_dot(self):
        if self.doc:
            return self.doc.split('\n')[0][-1] == '.'

    @property
    def deprecated(self):
        pattern = re.compile('.. deprecated:: ')
        return (self.method_name.startswith('pandas.Panel') or
                bool(pattern.search(self.summary)) or
                bool(pattern.search(self.extended_summary)))

    @property
    def mentioned_private_classes(self):
        return [klass for klass in PRIVATE_CLASSES if klass in self.raw_doc]

    @property
    def examples_errors(self):
        flags = doctest.NORMALIZE_WHITESPACE | doctest.IGNORE_EXCEPTION_DETAIL
        finder = doctest.DocTestFinder()
        runner = doctest.DocTestRunner(optionflags=flags)
        context = {'np': numpy, 'pd': pandas}
        error_msgs = ''
        for test in finder.find(self.raw_doc, self.method_name, globs=context):
            f = StringIO()
            runner.run(test, out=f.write)
            error_msgs += f.getvalue()
        return error_msgs


def get_api_items():
    api_fname = os.path.join(BASE_PATH, 'doc', 'source', 'api.rst')

    previous_line = current_section = current_subsection = ''
    position = None
    with open(api_fname) as f:
        for line in f:
            line = line.strip()
            if len(line) == len(previous_line):
                if set(line) == set('-'):
                    current_section = previous_line
                    continue
                if set(line) == set('~'):
                    current_subsection = previous_line
                    continue

            if line.startswith('.. currentmodule::'):
                current_module = line.replace('.. currentmodule::', '').strip()
                continue

            if line == '.. autosummary::':
                position = 'autosummary'
                continue

            if position == 'autosummary':
                if line == '':
                    position = 'items'
                    continue

            if position == 'items':
                if line == '':
                    position = None
                    continue
                item = line.strip()
                func = importlib.import_module(current_module)
                for part in item.split('.'):
                    func = getattr(func, part)

                yield ('.'.join([current_module, item]), func,
                       current_section, current_subsection)

            previous_line = line


def _csv_row(func_name, func_obj, section, subsection, in_api, seen={}):
    obj_type = type(func_obj).__name__
    original_callable = _to_original_callable(func_obj)
    if original_callable is None:
        return [func_name, obj_type] + [''] * 12, ''
    else:
        doc = Docstring(func_name, original_callable)
        key = doc.source_file_name, doc.source_file_def_line
        shared_code = seen.get(key, '')
        return [func_name,
                obj_type,
                in_api,
                int(doc.deprecated),
                section,
                subsection,
                doc.source_file_name,
                doc.source_file_def_line,
                doc.github_url,
                int(bool(doc.summary)),
                int(bool(doc.extended_summary)),
                int(doc.correct_parameters),
                int(bool(doc.examples)),
                shared_code], key


def validate_all():
    writer = csv.writer(sys.stdout)
    cols = ('Function or method',
            'Type',
            'In API doc',
            'Is deprecated',
            'Section',
            'Subsection',
            'File',
            'Code line',
            'GitHub link',
            'Has summary',
            'Has extended summary',
            'Parameters ok',
            'Has examples',
            'Shared code with')
    writer.writerow(cols)
    seen = {}
    api_items = list(get_api_items())
    for func_name, func, section, subsection in api_items:
        row, key = _csv_row(func_name, func, section, subsection,
                            in_api=1, seen=seen)
        seen[key] = func_name
        writer.writerow(row)

    api_item_names = set(list(zip(*api_items))[0])
    for class_ in (pandas.Series, pandas.DataFrame, pandas.Panel):
        for member in inspect.getmembers(class_):
            func_name = 'pandas.{}.{}'.format(class_.__name__, member[0])
            if (not member[0].startswith('_') and
                    func_name not in api_item_names):
                func = _load_obj(func_name)
                row, key = _csv_row(func_name, func, section='', subsection='',
                                    in_api=0)
                writer.writerow(row)

    return 0


def validate_one(func_name):
    """
    Validate the docstring for the given func_name

    Parameters
    ----------
    func_name : function
        Function whose docstring will be evaluated

    Returns
    -------
    int
        The number of errors found in the `func_name` docstring
    """
    func_obj = _load_obj(func_name)
    doc = Docstring(func_name, func_obj)

    sys.stderr.write(_output_header('Docstring ({})'.format(func_name)))
    sys.stderr.write('{}\n'.format(doc.clean_doc))

    errs = []
    wrns = []
    if doc.start_blank_lines != 1:
        errs.append('Docstring text (summary) should start in the line '
                    'immediately after the opening quotes (not in the same '
                    'line, or leaving a blank line in between)')
    if doc.end_blank_lines != 1:
        errs.append('Closing quotes should be placed in the line after '
                    'the last text in the docstring (do not close the '
                    'quotes in the same line as the text, or leave a '
                    'blank line between the last text and the quotes)')
    if doc.double_blank_lines:
        errs.append('Use only one blank line to separate sections or '
                    'paragraphs')

    if not doc.summary:
        errs.append('No summary found (a short summary in a single line '
                    'should be present at the beginning of the docstring)')
    else:
        if not doc.summary[0].isupper():
            errs.append('Summary does not start with a capital letter')
        if doc.summary[-1] != '.':
            errs.append('Summary does not end with a period')
        if (doc.is_function_or_method and
                doc.summary.split(' ')[0][-1] == 's'):
            errs.append('Summary must start with infinitive verb, '
                        'not third person (e.g. use "Generate" instead of '
                        '"Generates")')
    if not doc.extended_summary:
        wrns.append('No extended summary found')

    param_errs = doc.parameter_mismatches
    for param in doc.doc_parameters:
        if not param.startswith("*"):  # Check can ignore var / kwargs
            if not doc.parameter_type(param):
                param_errs.append('Parameter "{}" has no type'.format(param))
            else:
                if doc.parameter_type(param)[-1] == '.':
                    param_errs.append('Parameter "{}" type should '
                                      'not finish with "."'.format(param))

        if not doc.parameter_desc(param):
            param_errs.append('Parameter "{}" '
                              'has no description'.format(param))
        else:
            if not doc.parameter_desc(param)[0].isupper():
                param_errs.append('Parameter "{}" description '
                                  'should start with a '
                                  'capital letter'.format(param))
            if doc.parameter_desc(param)[-1] != '.':
                param_errs.append('Parameter "{}" description '
                                  'should finish with "."'.format(param))
    if param_errs:
        errs.append('Errors in parameters section')
        for param_err in param_errs:
            errs.append('\t{}'.format(param_err))

    if doc.is_function_or_method:
        if not doc.returns and "return" in doc.method_source:
            errs.append('No Returns section found')
        if not doc.yields and "yield" in doc.method_source:
            errs.append('No Yields section found')

    mentioned_errs = doc.mentioned_private_classes
    if mentioned_errs:
        errs.append('Private classes ({}) should not be mentioned in public '
                    'docstring.'.format(mentioned_errs))

    if not doc.see_also:
        wrns.append('See Also section not found')
    else:
        for rel_name, rel_desc in doc.see_also.items():
            if not rel_desc:
                errs.append('Missing description for '
                            'See Also "{}" reference'.format(rel_name))

    for line in doc.raw_doc.splitlines():
        if re.match("^ *\t", line):
            errs.append('Tabs found at the start of line "{}", '
                        'please use whitespace only'.format(line.lstrip()))

    examples_errs = ''
    if not doc.examples:
        wrns.append('No examples section found')
    else:
        examples_errs = doc.examples_errors
        if examples_errs:
            errs.append('Examples do not pass tests')

    sys.stderr.write(_output_header('Validation'))
    if errs:
        sys.stderr.write('Errors found:\n')
        for err in errs:
            sys.stderr.write('\t{}\n'.format(err))
    if wrns:
        sys.stderr.write('Warnings found:\n')
        for wrn in wrns:
            sys.stderr.write('\t{}\n'.format(wrn))

    if not errs:
        sys.stderr.write('Docstring for "{}" correct. :)\n'.format(func_name))

    if examples_errs:
        sys.stderr.write(_output_header('Doctests'))
        sys.stderr.write(examples_errs)

    return len(errs)


def main(function):
    if function is None:
        return validate_all()
    else:
        return validate_one(function)


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        description='validate pandas docstrings')
    argparser.add_argument('function',
                           nargs='?',
                           default=None,
                           help=('function or method to validate '
                                 '(e.g. pandas.DataFrame.head) '
                                 'if not provided, all docstrings '
                                 'are validated'))
    args = argparser.parse_args()
    sys.exit(main(args.function))
