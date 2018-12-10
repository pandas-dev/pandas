import argparse
import json
import os
import pickle
import re
import sys

import docutils.nodes

DOCTREE_PATH = '../doc/build/doctrees/{}.doctree'
RST_PATH = '../doc/source/{}.rst'

ERROR_MSGS = {
    'DT01': "Found bullet list within block quote",
    'WS01': "Indentation uses tabulation",
    'WS02': "Trailing whitespaces",
    'WS03': "Whitespaces in empty line",
}

STARTING_WHITESPACE_RE = re.compile(r'^(\s+).*\n$', re.MULTILINE)
TRAILING_WHITESPACE_RE = re.compile(r'^.*([\t ]+)\n$', re.MULTILINE)


class DocumentChecker(object):

    def __init__(self, raw_doc, doctree):
        self.doctree = doctree
        self.raw_doc = raw_doc

        self.issues = None

    def issue(self, code, match=None, **kwargs):
        """
        Parameters
        ----------
        code : str
            Error code.
        **kwargs
            Values for the variables in the error messages
        """
        issue = self.issues.setdefault(code, [])
        issue.append((self.find_line(match), kwargs))

    def find_line(self, match):
        if not match:
            return None
        lines = self.raw_doc[:match.start(0)].splitlines()
        return len(lines) + 1, match.group(0)

    def check_bullet_list_in_block_quote(self):
        for node in self.doctree.traverse(docutils.nodes.block_quote):
            match = node.first_child_matching_class(docutils.nodes.bullet_list)
            if match is not None:
                self.issue('DT01')

    def check_tabulator_as_indentation(self):
        matches = STARTING_WHITESPACE_RE.finditer(self.raw_doc)
        for match in matches:
            if '\t' in match.group(1):
                self.issue('WS01', match)

    def check_line_ends_with_whitespace(self):
        matches = TRAILING_WHITESPACE_RE.finditer(self.raw_doc)
        for match in matches:
            self.issue('WS02', match)

    def validate(self):
        self.issues = {}
        for func in dir(self):
            if func.startswith('check'):
                self.__class__.__dict__[func](self)

        return self.issues


def report(reports, output_format='default', errors=None):
    exit_status = 0
    if output_format == 'json':
        output = json.dumps(reports)
    else:
        if output_format == 'default':
            output_format = '{text}\n'
        elif output_format == 'azure':
            output_format = ('##vso[task.logissue type=error;'
                             'sourcepath={path};'
                             'linenumber={row};'
                             'code={code};'
                             ']{text}\n')
        else:
            raise ValueError('Unknown output_format "{}"'.format(
                output_format))

        output = ''
        for name, res in reports.items():
            for err_code, issues in res.items():
                # The script would be faster if instead of filtering the
                # errors after validating them, it didn't validate them
                # initially. But that would complicate the code too much
                if errors and err_code not in errors:
                    continue
                for issue, kwargs in issues:
                    exit_status += 1
                    row = issue[0] if issue else None
                    output += output_format.format(
                        name=name,
                        path=RST_PATH.format(name),
                        row=row,
                        code=err_code,
                        text='{}{}:: {}'.format(name,
                                                ':' + row if row else '',
                                                ERROR_MSGS[err_code]
                                                .format(kwargs)))

    sys.stdout.write(output)
    return exit_status


def validate(*docnames):
    result = {}
    for docname in docnames:
        with open(DOCTREE_PATH.format(docname), 'r+b') as file:
            doctree = pickle.load(file)

        with open(RST_PATH.format(docname), 'r') as file:
            raw_doc = file.read()

        checker = DocumentChecker(raw_doc, doctree)
        result[docname] = checker.validate()

    return result


def main(document, errors, output_format):
    if document:
        reports = validate(document)
    else:
        documentation_source = os.path.join(os.curdir, '../doc/source')
        docnames = []
        for root, dirs, files in os.walk(documentation_source):
            _, base_dir = root.split('../doc/source')
            for file in files:
                docname, ext = os.path.splitext(file)
                if not ext == '.rst':
                    continue
                docnames.append(os.path.join(base_dir, docname))
        reports = validate(*docnames)
    return report(reports, output_format=output_format, errors=errors)


if __name__ == '__main__':
    format_opts = 'default', 'json', 'azure'
    func_help = ('document to validate (e.g. io) '
                 'if not provided, all documents are validated')
    argparser = argparse.ArgumentParser(
        description='validate pandas documentation')
    add = argparser.add_argument
    add('document',
        nargs='?',
        default=None,
        help=func_help)
    add('--format', default='default', choices=format_opts,
        help='format of the output when validating '
             'multiple documents (ignored when validating one).'
             'It can be {}'.format(str(format_opts)[1:-1]))
    add('--errors', default=None,
        help='comma separated '
             'list of error codes to validate. By default it '
             'validates all errors (ignored when validating '
             'a single document)')

    args = argparser.parse_args()
    sys.exit(main(args.document,
                  args.errors.split(',') if args.errors else None,
                  args.format, ))
