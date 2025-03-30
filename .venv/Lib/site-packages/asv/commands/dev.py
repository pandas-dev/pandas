# Licensed under a 3-clause BSD style license - see LICENSE.rst

from .run import Run


class Dev(Run):
    @classmethod
    def setup_arguments(cls, subparsers):
        parser = subparsers.add_parser(
            "dev", help="Do a test run of a benchmark suite during development",
            description="""
                This runs a benchmark suite in a mode that is useful
                during development.  It is equivalent to
                ``asv run --python=same``""")

        cls._setup_arguments(parser, env_default_same=True)

        parser.set_defaults(func=cls.run_from_args)

        return parser

    @classmethod
    def run(cls, conf, **kwargs):
        raise RuntimeError(
            "`asv dev` has been removed. It was a shortcut of `asv run`, which "
            "you can use instead. It was removed because it caused confusion, "
            "in particular after changing what it was running in asv 0.5.x\n\n"
            "You can be interested in the next `asv run` arguments:\n"
            " - `--python=same`: to use the current environment and not create one\n"
            " - `--quick`: to run benchmarks just once\n"
            " - `--dry-run`: to not save the results\n"
        )
