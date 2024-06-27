# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
from os.path import dirname, join

from asv import config, commands, environment, util
from asv.commands import Command


def test_config():
    conf = config.Config.load(join(dirname(__file__), 'asv.conf.json'))

    assert conf.project == 'astropy'
    assert conf.matrix == {
        "numpy": ["1.8"],
        "Cython": [],
        "jinja2": [],
    }
    assert conf.benchmark_dir == 'benchmark'
    assert conf.branches == [None]
    assert conf.install_timeout == 3142  # gh-973
    assert conf.default_benchmark_timeout == 32.5  # gh-973


def test_config_default_benchmark_timeout():
    # gh-973
    conf = config.Config()
    assert conf.default_benchmark_timeout == 60.0


def test_config_default_install_timeout():
    # gh-391
    conf = config.Config()
    assert conf.install_timeout == 600


class CustomCommand(Command):
    @classmethod
    def setup_arguments(cls, subparsers):
        parser = subparsers.add_parser(
            "custom", help="Custom command",
            description="Juts a test.")

        parser.set_defaults(func=cls.run_from_args)

        return parser

    @classmethod
    def run_from_conf_args(cls, conf, args):
        pass


def test_load_plugin():
    os.chdir(dirname(__file__))

    parser, subparsers = commands.make_argparser()
    args = parser.parse_args(['custom'])

    assert hasattr(args, 'func')

    args.func(args)

    for env in util.iter_subclasses(environment.Environment):
        print(env.__name__)
        if env.__name__ == 'MyEnvironment':
            break
    else:
        assert False, "Custom plugin not loaded"
