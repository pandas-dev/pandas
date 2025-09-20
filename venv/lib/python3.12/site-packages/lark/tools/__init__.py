import sys
from argparse import ArgumentParser, FileType
from textwrap import indent
from logging import DEBUG, INFO, WARN, ERROR
from typing import Optional
import warnings

from lark import Lark, logger
try:
    from interegular import logger as interegular_logger
    has_interegular = True
except ImportError:
    has_interegular = False

lalr_argparser = ArgumentParser(add_help=False, epilog='Look at the Lark documentation for more info on the options')

flags = [
    ('d', 'debug'),
    'keep_all_tokens',
    'regex',
    'propagate_positions',
    'maybe_placeholders',
    'use_bytes'
]

options = ['start', 'lexer']

lalr_argparser.add_argument('-v', '--verbose', action='count', default=0, help="Increase Logger output level, up to three times")
lalr_argparser.add_argument('-s', '--start', action='append', default=[])
lalr_argparser.add_argument('-l', '--lexer', default='contextual', choices=('basic', 'contextual'))
lalr_argparser.add_argument('-o', '--out', type=FileType('w', encoding='utf-8'), default=sys.stdout, help='the output file (default=stdout)')
lalr_argparser.add_argument('grammar_file', type=FileType('r', encoding='utf-8'), help='A valid .lark file')

for flag in flags:
    if isinstance(flag, tuple):
        options.append(flag[1])
        lalr_argparser.add_argument('-' + flag[0], '--' + flag[1], action='store_true')
    elif isinstance(flag, str):
        options.append(flag)
        lalr_argparser.add_argument('--' + flag, action='store_true')
    else:
        raise NotImplementedError("flags must only contain strings or tuples of strings")


def build_lalr(namespace):
    logger.setLevel((ERROR, WARN, INFO, DEBUG)[min(namespace.verbose, 3)])
    if has_interegular:
        interegular_logger.setLevel(logger.getEffectiveLevel())
    if len(namespace.start) == 0:
        namespace.start.append('start')
    kwargs = {n: getattr(namespace, n) for n in options}
    return Lark(namespace.grammar_file, parser='lalr', **kwargs), namespace.out


def showwarning_as_comment(message, category, filename, lineno, file=None, line=None):
    # Based on warnings._showwarnmsg_impl
    text = warnings.formatwarning(message, category, filename, lineno, line)
    text = indent(text, '# ')
    if file is None:
        file = sys.stderr
        if file is None:
            return
    try:
        file.write(text)
    except OSError:
        pass


def make_warnings_comments():
    warnings.showwarning = showwarning_as_comment
