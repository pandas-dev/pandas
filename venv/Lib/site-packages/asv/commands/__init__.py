# Licensed under a 3-clause BSD style license - see LICENSE.rst

import argparse
import sys
import abc

from asv import config, util

from . import common_args

# This list is ordered in order of average workflow
command_order = [
    'Quickstart',
    'Machine',
    'Setup',
    'Run',
    'Dev',
    'Continuous',
    'Find',
    'Rm',
    'Publish',
    'Preview',
    'Profile',
    'Update',
    'Show',
    'Compare'
]


class Command(abc.ABC):
    @classmethod
    @abc.abstractmethod
    def setup_arguments(cls, subparsers):
        """
        Method to add the parser and arguments of the command.

        In most cases, `subparser.add_parser(cmd_name)` should be called to
        add the command (where `cmd_name` is for example `run` or `continuous`).
        And then call `add_argument` to the returned parser, to add the command
        arguments.
        """

    @classmethod
    def run_from_args(cls, args):
        from ..plugin_manager import plugin_manager
        conf = config.Config.load(args.config)
        for plugin in conf.plugins:
            plugin_manager.import_plugin(plugin)
        return cls.run_from_conf_args(conf, args)

    @classmethod
    @abc.abstractmethod
    def run_from_conf_args(cls, conf, args):
        """
        Call the run method (i.e. `cls.run`) with the right arguments from
        `conf` and `args`. In most cases this is just something like
        `clf.run(conf=conf, my_arg_1=args.my_arg_1, my_arg_2=args.my_arg_2, ...)`.
        """


def make_argparser():
    """
    The top-level entry point for the asv script.

    Most of the real work is handled by the subcommands in the
    commands subpackage.
    """
    def help(args):
        parser.print_help()
        sys.exit(0)

    parser = argparse.ArgumentParser(
        "asv",
        description="Airspeed Velocity: Simple benchmarking tool for Python")

    common_args.add_global_arguments(parser, suppress_defaults=False)

    subparsers = parser.add_subparsers(
        title='subcommands',
        description='valid subcommands')

    help_parser = subparsers.add_parser(
        "help", help="Display usage information")
    help_parser.set_defaults(func=help)

    commands = dict((x.__name__, x) for x in util.iter_subclasses(Command))

    for command in command_order:
        subparser = commands[str(command)].setup_arguments(subparsers)
        common_args.add_global_arguments(subparser)
        del commands[command]

    for name, command in sorted(commands.items()):
        subparser = command.setup_arguments(subparsers)
        common_args.add_global_arguments(subparser)

    return parser, subparsers


def _make_docstring():
    parser, subparsers = make_argparser()

    lines = []

    for p in subparsers.choices.values():
        lines.append(f".. _cmd-{p.prog.replace(' ', '-')}:")
        lines.append('')
        lines.append(p.prog)
        lines.append('-' * len(p.prog))
        lines.append('::')
        lines.append('')
        lines.extend('   ' + x for x in p.format_help().splitlines())
        lines.append('')

    return '\n'.join(lines)
