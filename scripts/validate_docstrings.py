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
import json
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


def get_api_items():
    """
    Parse api.rst file from the documentation, and extract all the functions,
    methods, classes, attributes... This should include all pandas public API.

    Yields
    ------
    name : str
        The name of the object (e.g. 'pandas.Series.str.upper).
    func : function
        The object itself. In most cases this will be a function or method,
        but it can also be classes, properties, cython objects...
    section : str
        The name of the section in the API page where the object item is
        located.
    subsection : str
        The name of the subsection in the API page where the object item is
        located.
    """
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


class Docstring(object):
    def __init__(self, name):
        self.name = name
        obj = self._load_obj(name)
        self.obj = obj
        self.code_obj = self._to_original_callable(obj)
        self.raw_doc = obj.__doc__ or ''
        self.clean_doc = pydoc.getdoc(obj)
        self.doc = NumpyDocString(self.clean_doc)

    def __len__(self):
        return len(self.raw_doc)

    @staticmethod
    def _load_obj(name):
        """
        Import Python object from its name as string.

        Parameters
        ----------
        name : str
            Object name to import (e.g. pandas.Series.str.upper)

        Returns
        -------
        object
            Python object that can be a class, method, function...

        Examples
        --------
        >>> Docstring._load_obj('pandas.Series')
        <class 'pandas.core.series.Series'>
        """
        for maxsplit in range(1, name.count('.') + 1):
            # TODO when py3 only replace by: module, *func_parts = ...
            func_name_split = name.rsplit('.', maxsplit)
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
                              'from "{}"'.format(name))

        for part in func_parts:
            obj = getattr(obj, part)
        return obj

    @staticmethod
    def _to_original_callable(obj):
        """
        Find the Python object that contains the source code ot the object.

        This is useful to find the place in the source code (file and line
        number) where a docstring is defined. It does not currently work for
        all cases, but it should help find some (properties...).
        """
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

    @property
    def type(self):
        return type(self.obj).__name__

    @property
    def is_function_or_method(self):
        # TODO(py27): remove ismethod
        return (inspect.isfunction(self.obj)
                or inspect.ismethod(self.obj))

    @property
    def source_file_name(self):
        """
        File name where the object is implemented (e.g. pandas/core/frame.py).
        """
        try:
            fname = inspect.getsourcefile(self.code_obj)
        except TypeError:
            # In some cases the object is something complex like a cython
            # object that can't be easily introspected. An it's better to
            # return the source code file of the object as None, than crash
            pass
        else:
            if fname:
                fname = os.path.relpath(fname, BASE_PATH)
                return fname

    @property
    def source_file_def_line(self):
        """
        Number of line where the object is defined in its file.
        """
        try:
            return inspect.getsourcelines(self.code_obj)[-1]
        except (OSError, TypeError):
            # In some cases the object is something complex like a cython
            # object that can't be easily introspected. An it's better to
            # return the line number as None, than crash
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
        if inspect.isclass(self.obj):
            if hasattr(self.obj, '_accessors') and (
                    self.name.split('.')[-1] in
                    self.obj._accessors):
                # accessor classes have a signature but don't want to show this
                return tuple()
        try:
            sig = signature(self.obj)
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
        try:
            return inspect.getsource(self.obj)
        except TypeError:
            return ''

    @property
    def first_line_ends_in_dot(self):
        if self.doc:
            return self.doc.split('\n')[0][-1] == '.'

    @property
    def deprecated(self):
        pattern = re.compile('.. deprecated:: ')
        return (self.name.startswith('pandas.Panel') or
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
        for test in finder.find(self.raw_doc, self.name, globs=context):
            f = StringIO()
            runner.run(test, out=f.write)
            error_msgs += f.getvalue()
        return error_msgs


def validate_one(func_name):
    """
    Validate the docstring for the given func_name

    Parameters
    ----------
    func_name : function
        Function whose docstring will be evaluated (e.g. pandas.read_csv).

    Returns
    -------
    dict
        A dictionary containing all the information obtained from validating
        the docstring.
    """
    doc = Docstring(func_name)

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
        if doc.summary != doc.summary.lstrip():
            errs.append('Summary contains heading whitespaces.')
        elif (doc.is_function_or_method and
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

    return {'type': doc.type,
            'docstring': doc.clean_doc,
            'deprecated': doc.deprecated,
            'file': doc.source_file_name,
            'file_line': doc.source_file_def_line,
            'github_link': doc.github_url,
            'errors': errs,
            'warnings': wrns,
            'examples_errors': examples_errs}


def validate_all():
    """
    Execute the validation of all docstrings, and return a dict with the
    results.

    Returns
    -------
    dict
        A dictionary with an item for every function/method... containing
        all the validation information.
    """
    result = {}
    seen = {}

    # functions from the API docs
    api_items = list(get_api_items())
    for func_name, func_obj, section, subsection in api_items:
        doc_info = validate_one(func_name)
        result[func_name] = doc_info

        shared_code_key = doc_info['file'], doc_info['file_line']
        shared_code = seen.get(shared_code_key, '')
        result[func_name].update({'in_api': True,
                                  'section': section,
                                  'subsection': subsection,
                                  'shared_code_with': shared_code})

        seen[shared_code_key] = func_name

    # functions from introspecting Series, DataFrame and Panel
    api_item_names = set(list(zip(*api_items))[0])
    for class_ in (pandas.Series, pandas.DataFrame, pandas.Panel):
        for member in inspect.getmembers(class_):
            func_name = 'pandas.{}.{}'.format(class_.__name__, member[0])
            if (not member[0].startswith('_') and
                    func_name not in api_item_names):
                doc_info = validate_one(func_name)
                result[func_name] = doc_info
                result[func_name]['in_api'] = False

    return result


def main(func_name, fd):
    def header(title, width=80, char='#'):
        full_line = char * width
        side_len = (width - len(title) - 2) // 2
        adj = '' if len(title) % 2 == 0 else ' '
        title_line = '{side} {title}{adj} {side}'.format(side=char * side_len,
                                                         title=title,
                                                         adj=adj)

        return '\n{full_line}\n{title_line}\n{full_line}\n\n'.format(
            full_line=full_line, title_line=title_line)

    if func_name is None:
        json_doc = validate_all()
        fd.write(json.dumps(json_doc))
    else:
        doc_info = validate_one(func_name)

        fd.write(header('Docstring ({})'.format(func_name)))
        fd.write('{}\n'.format(doc_info['docstring']))
        fd.write(header('Validation'))
        if doc_info['errors']:
            fd.write('Errors found:\n')
            for err in doc_info['errors']:
                fd.write('\t{}\n'.format(err))
        if doc_info['warnings']:
            fd.write('Warnings found:\n')
            for wrn in doc_info['warnings']:
                fd.write('\t{}\n'.format(wrn))

        if not doc_info['errors']:
            fd.write('Docstring for "{}" correct. :)\n'.format(func_name))

        if doc_info['examples_errors']:
            fd.write(header('Doctests'))
            fd.write(doc_info['examples_errors'])


if __name__ == '__main__':
    func_help = ('function or method to validate (e.g. pandas.DataFrame.head) '
                 'if not provided, all docstrings are validated and returned '
                 'as JSON')
    argparser = argparse.ArgumentParser(
        description='validate pandas docstrings')
    argparser.add_argument('function',
                           nargs='?',
                           default=None,
                           help=func_help)
    args = argparser.parse_args()
    sys.exit(main(args.function, sys.stdout))
