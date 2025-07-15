
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from ..benchmarks import Benchmarks
from . import Command, common_args
from ..repo import get_repo
from .. import environment


class Check(Command):
    @classmethod
    def setup_arguments(cls, subparsers):
        parser = subparsers.add_parser(
            "check", help="Import and check benchmark suite, but do not run benchmarks",
            description="""
                This imports and checks basic validity of the benchmark suite, but
                does not run the benchmark target code""")

        common_args.add_environment(parser, default_same=False)
        parser.set_defaults(func=cls.run_from_args)

        return parser

    @classmethod
    def run_from_conf_args(cls, conf, args, **kwargs):
        return cls.run(conf, env_spec=args.env_spec, **kwargs)

    @classmethod
    def run(cls, conf, env_spec=None):
        environments = list(environment.get_environments(conf, env_spec))

        if environment.is_existing_only(environments):
            # No repository required, so skip using it
            conf.dvcs = "none"

        repo = get_repo(conf)

        commit_hashes = []  # use fallback defaults
        Benchmarks.discover(conf, repo, environments, commit_hashes, check=True)
        return 0
