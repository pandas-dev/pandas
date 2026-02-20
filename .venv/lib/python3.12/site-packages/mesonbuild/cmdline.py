# SPDX-License-Identifier: Apache-2.0
# Copyright 2013-2025 The Meson development team

from __future__ import annotations

import argparse
import ast
import configparser
import os
import shlex
import typing as T
from itertools import chain

from . import options
from .mesonlib import MesonException
from .options import OptionKey

if T.TYPE_CHECKING:
    from typing_extensions import Protocol

    # typeshed
    StrOrBytesPath = T.Union[str, bytes, os.PathLike[str], os.PathLike[bytes]]

    class SharedCMDOptions(Protocol):
        """Representation of command line options from Meson setup, configure,
        and dist.

        :param cmd_line_options: command line options parsed into an OptionKey:
            str mapping
        """

        cmd_line_options: T.Dict[OptionKey, T.Optional[str]]
        cross_file: T.List[str]
        native_file: T.List[str]


class CmdLineFileParser(configparser.ConfigParser):
    def __init__(self) -> None:
        # We don't want ':' as key delimiter, otherwise it would break when
        # storing subproject options like "subproject:option=value"
        super().__init__(delimiters=['='], interpolation=None)

    def read(self, filenames: T.Union['StrOrBytesPath', T.Iterable['StrOrBytesPath']], encoding: T.Optional[str] = 'utf-8') -> T.List[str]:
        return super().read(filenames, encoding)

    def optionxform(self, optionstr: str) -> str:
        # Don't call str.lower() on keys
        return optionstr


def get_cmd_line_file(build_dir: str) -> str:
    return os.path.join(build_dir, 'meson-private', 'cmd_line.txt')

def read_cmd_line_file(build_dir: str, options: SharedCMDOptions) -> None:
    filename = get_cmd_line_file(build_dir)
    if not os.path.isfile(filename):
        return

    config = CmdLineFileParser()
    config.read(filename)

    # Do a copy because config is not really a dict. options.cmd_line_options
    # overrides values from the file.
    d = {OptionKey.from_string(k): v for k, v in config['options'].items()}
    d.update(options.cmd_line_options)
    options.cmd_line_options = d

    properties = config['properties']
    if not options.cross_file:
        options.cross_file = ast.literal_eval(properties.get('cross_file', '[]'))
    if not options.native_file:
        # This will be a string in the form: "['first', 'second', ...]", use
        # literal_eval to get it into the list of strings.
        options.native_file = ast.literal_eval(properties.get('native_file', '[]'))

def write_cmd_line_file(build_dir: str, options: SharedCMDOptions) -> None:
    filename = get_cmd_line_file(build_dir)
    config = CmdLineFileParser()

    properties: T.Dict[str, T.List[str]] = {}
    if options.cross_file:
        properties['cross_file'] = options.cross_file
    if options.native_file:
        properties['native_file'] = options.native_file

    config['options'] = {str(k): str(v) for k, v in options.cmd_line_options.items()}
    config['properties'] = {k: repr(v) for k, v in properties.items()}
    with open(filename, 'w', encoding='utf-8') as f:
        config.write(f)

def update_cmd_line_file(build_dir: str, options: SharedCMDOptions) -> None:
    filename = get_cmd_line_file(build_dir)
    config = CmdLineFileParser()
    config.read(filename)
    for k, v in options.cmd_line_options.items():
        keystr = str(k)
        if v is not None:
            config['options'][keystr] = str(v)
        elif keystr in config['options']:
            del config['options'][keystr]

    with open(filename, 'w', encoding='utf-8') as f:
        config.write(f)

def format_cmd_line_options(options: SharedCMDOptions) -> str:
    cmdline = ['-D{}={}'.format(str(k), v) for k, v in options.cmd_line_options.items()]
    if options.cross_file:
        cmdline += [f'--cross-file={f}' for f in options.cross_file]
    if options.native_file:
        cmdline += [f'--native-file={f}' for f in options.native_file]
    return ' '.join([shlex.quote(x) for x in cmdline])


class KeyNoneAction(argparse.Action):
    """
    Custom argparse Action that stores values in a dictionary as keys with value None.
    """

    def __init__(self, option_strings: str, dest: str, nargs: T.Optional[T.Union[int, str]] = None, **kwargs: T.Any) -> None:
        assert nargs is None or nargs == 1
        super().__init__(option_strings, dest, nargs=1, **kwargs)

    def __call__(self, parser: argparse.ArgumentParser, namespace: argparse.Namespace,
                 arg: T.List[str], option_string: str = None) -> None: # type: ignore[override]
        current_dict = getattr(namespace, self.dest)
        if current_dict is None:
            current_dict = {}
            setattr(namespace, self.dest, current_dict)

        key = OptionKey.from_string(arg[0])
        current_dict[key] = None


class KeyValueAction(argparse.Action):
    """
    Custom argparse Action that parses KEY=VAL arguments and stores them in a dictionary.
    """

    def __init__(self, option_strings: str, dest: str, nargs: T.Optional[T.Union[int, str]] = None, **kwargs: T.Any) -> None:
        assert nargs is None or nargs == 1
        super().__init__(option_strings, dest, nargs=1, **kwargs)

    def __call__(self, parser: argparse.ArgumentParser, namespace: argparse.Namespace,
                 arg: T.List[str], option_string: str = None) -> None: # type: ignore[override]
        current_dict = getattr(namespace, self.dest)
        if current_dict is None:
            current_dict = {}
            setattr(namespace, self.dest, current_dict)

        try:
            keystr, value = arg[0].split('=', 1)
            key = OptionKey.from_string(keystr)
            current_dict[key] = value
        except ValueError:
            parser.error(f'The argument for option {option_string!r} must be in OPTION=VALUE format.')


def register_builtin_arguments(parser: argparse.ArgumentParser) -> None:
    for n, b in options.BUILTIN_OPTIONS.items():
        options.option_to_argparse(b, n, parser, '')
    for n, b in options.BUILTIN_OPTIONS_PER_MACHINE.items():
        options.option_to_argparse(b, n, parser, ' (just for host machine)')
        options.option_to_argparse(b, n.as_build(), parser, ' (just for build machine)')
    parser.add_argument('-D', action=KeyValueAction, dest='cmd_line_options', default={}, metavar="option=value",
                        help='Set the value of an option, can be used several times to set multiple options.')

def parse_cmd_line_options(args: SharedCMDOptions) -> None:
    # Merge builtin options set with --option into the dict.
    for key in chain(
            options.BUILTIN_OPTIONS.keys(),
            (k.as_build() for k in options.BUILTIN_OPTIONS_PER_MACHINE.keys()),
            options.BUILTIN_OPTIONS_PER_MACHINE.keys(),
    ):
        name = str(key)
        value = getattr(args, name, None)
        if value is not None:
            if key in args.cmd_line_options:
                cmdline_name = options.argparse_name_to_arg(name)
                raise MesonException(
                    f'Got argument {name} as both -D{name} and {cmdline_name}. Pick one.')
            args.cmd_line_options[key] = value
            delattr(args, name)
