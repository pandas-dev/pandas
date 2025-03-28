# Licensed under a 3-clause BSD style license - see LICENSE.rst

import contextlib
import io
import os
import pstats
import sys
import tempfile

from asv_runner.console import color_print

from . import Command, common_args
from ..benchmarks import Benchmarks
from ..console import log
from ..environment import get_environments, is_existing_only
from ..machine import Machine
from ..profiling import ProfilerGui
from ..repo import get_repo, NoSuchNameError
from ..results import iter_results_for_machine
from ..runner import run_benchmarks
from ..util import hash_equal, iter_subclasses
from .. import util


@contextlib.contextmanager
def temp_profile(profile_data):
    profile_fd, profile_path = tempfile.mkstemp()
    try:
        with io.open(profile_fd, 'wb', closefd=True) as fd:
            fd.write(profile_data)

        yield profile_path
    finally:
        os.remove(profile_path)


class Profile(Command):
    @classmethod
    def setup_arguments(cls, subparsers):
        parser = subparsers.add_parser(
            "profile",
            help="""Run the profiler on a particular benchmark on a
            particular revision""",
            description="Profile a benchmark")

        parser.add_argument(
            'benchmark',
            help="""The benchmark to profile.  Must be a
            fully-specified benchmark name. For parameterized benchmark, it
            must include the parameter combination to use, e.g.:
            benchmark_name\\(param0, param1, ...\\)""")
        parser.add_argument(
            'revision', nargs='?',
            help="""The revision of the project to profile.  May be a
            commit hash, or a tag or branch name.""")
        parser.add_argument(
            '--gui', '-g',
            help="""Display the profile in the given gui.  Use
            --gui=list to list available guis.""")
        parser.add_argument(
            '--output', '-o',
            help="""Save the profiling information to the given file.
            This file is in the format written by the `cProfile`
            standard library module.  If not provided, prints a simple
            text-based profiling report to the console.""")
        parser.add_argument(
            '--force', '-f', action='store_true',
            help="""Forcibly re-run the profile, even if the data
            already exists in the results database.""")
        common_args.add_environment(parser)
        common_args.add_launch_method(parser)

        parser.set_defaults(func=cls.run_from_args)

        return parser

    @classmethod
    def find_guis(cls):
        cls.guis = {}
        for x in iter_subclasses(ProfilerGui):
            if x.name is not None and x.is_available():
                cls.guis[x.name] = x

    @classmethod
    def run_from_conf_args(cls, conf, args, **kwargs):
        return cls.run(
            conf=conf, benchmark=args.benchmark, revision=args.revision,
            gui=args.gui, output=args.output, force=args.force,
            env_spec=args.env_spec, launch_method=args.launch_method,
            **kwargs)

    @classmethod
    def run(cls, conf, benchmark, revision=None, gui=None, output=None,
            force=False, env_spec=None, launch_method=None,
            _machine_file=None):
        cls.find_guis()

        if gui == 'list':
            log.info("Available profiler GUIs:")
            with log.indent():
                for x in cls.guis.values():
                    log.info(f"{x.name}: {x.description}")
            return

        if gui is not None and gui not in cls.guis:
            raise util.UserError(
                f"Unknown profiler GUI {gui}")

        if benchmark is None:
            raise util.UserError(
                "Must specify benchmark to run")

        environments = list(get_environments(conf, env_spec))

        if is_existing_only(environments):
            # No repository required, so skip using it
            conf.dvcs = "none"

        repo = get_repo(conf)
        repo.pull()

        machine_name = Machine.get_unique_machine_name()
        if revision is None:
            rev = conf.branches[0]
        else:
            rev = revision

        try:
            commit_hash = repo.get_hash_from_name(rev)
        except NoSuchNameError as exc:
            raise util.UserError(f"Unknown commit {exc}")

        profile_data = None
        checked_out = set()

        # First, we see if we already have the profile in the results
        # database
        env = None
        if not force and commit_hash:
            for result in iter_results_for_machine(
                    conf.results_dir, machine_name):
                if hash_equal(commit_hash, result.commit_hash):
                    if result.has_profile(benchmark):
                        # Only take the first one
                        env_matched = util.get_matching_environment(environments, result)

                        if env_matched:
                            if result.env_name not in checked_out:
                                # We need to checkout the correct commit so that
                                # the line numbers in the profile data match up with
                                # what's in the source tree.
                                env_matched.checkout_project(repo, commit_hash)
                                checked_out.add(result.env_name)
                            profile_data = result.get_profile(benchmark)
                            env = env_matched
                            break

        if profile_data is None:
            if len(environments) == 0:
                log.error("No environments selected")
                return

            if revision is not None:
                for env in environments:
                    if not env.can_install_project():
                        raise util.UserError(
                            "An explicit revision may not be specified when "
                            "using an existing environment.")

            if env is None:
                # Fallback
                env = environments[0]

            if env.python != "{0}.{1}".format(*sys.version_info[:2]):
                raise util.UserError(
                    "Profiles must be run in the same version of Python as the "
                    "asv main process")

            benchmarks = Benchmarks.discover(conf, repo, environments,
                                             [commit_hash],
                                             regex=f'^{benchmark}$')

            if len(benchmarks) == 0:
                raise util.UserError(f"'{benchmark}' benchmark not found")
            elif len(benchmarks) > 1:
                exact_matches = benchmarks.filter_out([x for x in benchmarks if x != benchmark])
                if len(exact_matches) == 1:
                    log.warning("'{0}' matches more than one benchmark, "
                                "using exact match".format(benchmark))
                    benchmarks = exact_matches
                else:
                    raise util.UserError(f"'{benchmark}' matches more than one benchmark")

            benchmark_name, = benchmarks.keys()

            if not force:
                log.info(
                    "Profile data does not already exist. "
                    "Running profiler now.")
            else:
                log.info("Running profiler")
            with log.indent():
                env.install_project(conf, repo, commit_hash)

                results = run_benchmarks(
                    benchmarks, env, show_stderr=True, quick=False, profile=True,
                    launch_method=launch_method)

                profile_data = results.get_profile(benchmark_name)

        log.flush()

        if gui is not None:
            log.debug(f"Opening gui {gui}")
            with temp_profile(profile_data) as profile_path:
                return cls.guis[gui].open_profiler_gui(profile_path)
        elif output is not None:
            with io.open(output, 'wb') as fd:
                fd.write(profile_data)
        else:
            color_print('')
            with temp_profile(profile_data) as profile_path:
                stats = pstats.Stats(profile_path)
                stats.strip_dirs()  # Addresses gh-71
                stats.sort_stats('cumulative')
                stats.print_stats()
