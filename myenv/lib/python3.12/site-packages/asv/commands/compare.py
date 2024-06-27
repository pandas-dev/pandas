# Licensed under a 3-clause BSD style license - see LICENSE.rst

import itertools
import math

import tabulate
from asv_runner.console import color_print
from asv_runner.statistics import get_err

from . import Command, common_args
from ..benchmarks import Benchmarks
from ..machine import iter_machine_files
from ..results import iter_results_for_machine_and_hash
from ..repo import get_repo, NoSuchNameError
from ..util import human_value, load_json
from ..console import log
from ..environment import get_environments
from .. import util, statistics


def mean(values):
    if all([value is None for value in values]):
        return None
    else:
        values = [value for value in values if value is not None]
        return sum(values) / float(len(values))


def unroll_result(benchmark_name, params, *values):
    """
    Iterate through parameterized result values

    Yields
    ------
    name
        Strings of the form "benchmark_name(value1, value2)" with
        parameter values substituted in. For non-parameterized
        results, simply the benchmark name.
    value
        Benchmark timing or other scalar value.

    """
    num_comb = 1
    for p in params:
        num_comb *= len(p)

    values = list(values)
    for j in range(len(values)):
        if values[j] is None:
            values[j] = [None] * num_comb

    for params, value in zip(itertools.product(*params), zip(*values)):
        if params == ():
            name = benchmark_name
        else:
            name = f"{benchmark_name}({', '.join(params)})"
        yield (name,) + value


def _isna(value):
    # None (failed) or NaN (skipped)
    return value is None or value != value


def _is_result_better(a, b, a_ss, b_ss, factor, use_stats=True):
    """
    Check if result 'a' is better than 'b' by the given factor,
    possibly taking confidence intervals into account.

    """

    if use_stats and a_ss and b_ss and a_ss[0] and b_ss[0] and (
            a_ss[0].get('repeat', 0) != 1 and b_ss[0].get('repeat', 0) != 1):
        # Return False if estimates don't differ.
        #
        # Special-case the situation with only one sample, in which
        # case we do the comparison only based on `factor` as there's
        # not enough data to do statistics.
        if not statistics.is_different(a_ss[1], b_ss[1],
                                       a_ss[0], b_ss[0]):
            return False

    return a < b / factor


class Compare(Command):

    @classmethod
    def setup_arguments(cls, subparsers):
        parser = subparsers.add_parser(
            "compare",
            help="""Compare the benchmark results between two revisions
                    (averaged over configurations)""",
            description="Compare two sets of results")

        parser.add_argument(
            'revision1',
            help="""The reference revision.""")

        parser.add_argument(
            'revision2',
            help="""The revision being compared.""")

        common_args.add_compare(parser, sort_default='default', only_changed_default=False)

        parser.add_argument(
            '--machine', '-m', type=str, default=None,
            help="""The machine to compare the revisions for.""")

        common_args.add_environment(parser)

        parser.set_defaults(func=cls.run_from_args)

        return parser

    @classmethod
    def run_from_conf_args(cls, conf, args):
        return cls.run(conf=conf,
                       hash_1=args.revision1,
                       hash_2=args.revision2,
                       factor=args.factor, split=args.split,
                       only_changed=args.only_changed, sort=args.sort,
                       machine=args.machine,
                       env_spec=args.env_spec,
                       use_stats=args.use_stats)

    @classmethod
    def run(cls, conf, hash_1, hash_2, factor=None, split=False, only_changed=False,
            sort='default', machine=None, env_spec=None, use_stats=True):

        repo = get_repo(conf)
        try:
            hash_1 = repo.get_hash_from_name(hash_1)
        except NoSuchNameError:
            pass

        try:
            hash_2 = repo.get_hash_from_name(hash_2)
        except NoSuchNameError:
            pass

        if env_spec:
            env_names = ([env.name for env in get_environments(conf, env_spec, verbose=False)] +
                         list(env_spec))
        else:
            env_names = None

        machines = []
        for path in iter_machine_files(conf.results_dir):
            d = load_json(path)
            machines.append(d['machine'])

        if len(machines) == 0:
            raise util.UserError("No results found")
        elif machine is None:
            if len(machines) > 1:
                raise util.UserError(
                    "Results available for several machines: {0} - "
                    "specify which one to use with the --machine option".format(
                        '/'.join(machines)))
            else:
                machine = machines[0]
        elif machine not in machines:
            raise util.UserError(
                f"Results for machine '{machine} not found")

        commit_names = {hash_1: repo.get_name_from_hash(hash_1),
                        hash_2: repo.get_name_from_hash(hash_2)}

        cls.print_table(conf, hash_1, hash_2, factor=factor, split=split, use_stats=use_stats,
                        only_changed=only_changed, sort=sort,
                        machine=machine, env_names=env_names, commit_names=commit_names)

    @classmethod
    def print_table(cls, conf, hash_1, hash_2, factor, split,
                    resultset_1=None, resultset_2=None, machine=None,
                    only_changed=False, sort='default',
                    use_stats=True, env_names=None, commit_names=None):
        results_1 = {}
        results_2 = {}
        ss_1 = {}
        ss_2 = {}
        versions_1 = {}
        versions_2 = {}
        units = {}

        benchmarks = Benchmarks.load(conf)

        if commit_names is None:
            commit_names = {}

        def results_default_iter(commit_hash):
            for result in iter_results_for_machine_and_hash(
                    conf.results_dir, machine, commit_hash):
                if env_names is not None and result.env_name not in env_names:
                    continue
                for key in result.get_all_result_keys():
                    params = result.get_result_params(key)
                    result_value = result.get_result_value(key, params)
                    result_stats = result.get_result_stats(key, params)
                    result_samples = result.get_result_samples(key, params)
                    result_version = result.benchmark_version.get(key)
                    yield (key, params, result_value, result_stats, result_samples,
                           result_version, result.params['machine'], result.env_name)

        if resultset_1 is None:
            resultset_1 = results_default_iter(hash_1)

        if resultset_2 is None:
            resultset_2 = results_default_iter(hash_2)

        machine_env_names = set()

        for key, params, value, stats, samples, version, machine, env_name in resultset_1:
            machine_env_name = f"{machine}/{env_name}"
            machine_env_names.add(machine_env_name)
            for name, value, stats, samples in unroll_result(key, params, value, stats, samples):
                units[(name, machine_env_name)] = benchmarks.get(key, {}).get('unit')
                results_1[(name, machine_env_name)] = value
                ss_1[(name, machine_env_name)] = (stats, samples)
                versions_1[(name, machine_env_name)] = version

        for key, params, value, stats, samples, version, machine, env_name in resultset_2:
            machine_env_name = f"{machine}/{env_name}"
            machine_env_names.add(machine_env_name)
            for name, value, stats, samples in unroll_result(key, params, value, stats, samples):
                units[(name, machine_env_name)] = benchmarks.get(key, {}).get('unit')
                results_2[(name, machine_env_name)] = value
                ss_2[(name, machine_env_name)] = (stats, samples)
                versions_2[(name, machine_env_name)] = version

        if len(results_1) == 0:
            raise util.UserError(
                f"Did not find results for commit {hash_1}")

        if len(results_2) == 0:
            raise util.UserError(
                f"Did not find results for commit {hash_2}")

        benchmarks_1 = set(results_1.keys())
        benchmarks_2 = set(results_2.keys())

        joint_benchmarks = sorted(list(benchmarks_1 | benchmarks_2))

        bench = {}

        if split:
            bench['green'] = []
            bench['red'] = []
            bench['lightgrey'] = []
            bench['default'] = []
        else:
            bench['all'] = []

        worsened = False
        improved = False

        for benchmark in joint_benchmarks:
            if benchmark in results_1:
                time_1 = results_1[benchmark]
            else:
                time_1 = math.nan

            if benchmark in results_2:
                time_2 = results_2[benchmark]
            else:
                time_2 = math.nan

            if benchmark in ss_1 and ss_1[benchmark][0]:
                err_1 = get_err(time_1, ss_1[benchmark][0])
            else:
                err_1 = None

            if benchmark in ss_2 and ss_2[benchmark][0]:
                err_2 = get_err(time_2, ss_2[benchmark][0])
            else:
                err_2 = None

            version_1 = versions_1.get(benchmark)
            version_2 = versions_2.get(benchmark)

            if _isna(time_1) or _isna(time_2):
                ratio = 'n/a'
                ratio_num = 1e9
            else:
                try:
                    ratio_num = time_2 / time_1
                    ratio = f"{ratio_num:6.2f}"
                except ZeroDivisionError:
                    ratio_num = 1e9
                    ratio = "n/a"

            if (version_1 is not None and version_2 is not None and
                    version_1 != version_2):
                # not comparable
                color = 'lightgrey'
                mark = 'x'
            elif time_1 is not None and time_2 is None:
                # introduced a failure
                color = 'red'
                mark = '!'
                worsened = True
            elif time_1 is None and time_2 is not None:
                # fixed a failure
                color = 'green'
                mark = ' '
                improved = True
            elif time_1 is None and time_2 is None:
                # both failed
                color = 'default'
                mark = ' '
            elif _isna(time_1) or _isna(time_2):
                # either one was skipped
                color = 'default'
                mark = ' '
            elif _is_result_better(time_2, time_1,
                                   ss_2.get(benchmark), ss_1.get(benchmark),
                                   factor, use_stats=use_stats):
                color = 'green'
                mark = '-'
                improved = True
            elif _is_result_better(time_1, time_2,
                                   ss_1.get(benchmark), ss_2.get(benchmark),
                                   factor, use_stats=use_stats):
                color = 'red'
                mark = '+'
                worsened = True
            else:
                color = 'default'
                mark = ' '

                # Mark statistically insignificant results
                if (_is_result_better(time_1, time_2, None, None, factor) or
                        _is_result_better(time_2, time_1, None, None, factor)):
                    ratio = "~" + ratio.strip()

            if only_changed and mark in (' ', 'x'):
                continue

            unit = units[benchmark]

            details = "{0:1s} {1:>15s}  {2:>15s} {3:>8s}  ".format(
                mark,
                human_value(time_1, unit, err=err_1),
                human_value(time_2, unit, err=err_2),
                ratio)
            split_line = details.split()
            if len(machine_env_names) > 1:
                benchmark_name = "{} [{}]".format(*benchmark)
            else:
                benchmark_name = benchmark[0]
            if len(split_line) == 4:
                split_line += [benchmark_name]
            else:
                split_line = [' '] + split_line + [benchmark_name]
            if split:
                bench[color].append(split_line)
            else:
                bench['all'].append(split_line)

        if split:
            keys = ['green', 'default', 'red', 'lightgrey']
        else:
            keys = ['all']

        titles = {}
        titles['green'] = "Benchmarks that have improved:"
        titles['default'] = "Benchmarks that have stayed the same:"
        titles['red'] = "Benchmarks that have got worse:"
        titles['lightgrey'] = "Benchmarks that are not comparable:"
        titles['all'] = "All benchmarks:"

        log.flush()

        for key in keys:

            if len(bench[key]) == 0:
                continue

            if not only_changed:
                color_print("")
                color_print(titles[key])
                color_print("")

            name_1 = commit_names.get(hash_1)
            if name_1:
                name_1 = f'<{name_1}>'
            else:
                name_1 = ''

            name_2 = commit_names.get(hash_2)
            if name_2:
                name_2 = f'<{name_2}>'
            else:
                name_2 = ''

            if sort == 'default':
                pass
            elif sort == 'ratio':
                bench[key].sort(key=lambda v: v[3], reverse=True)
            elif sort == 'name':
                bench[key].sort(key=lambda v: v[2])
            else:
                raise ValueError("Unknown 'sort'")

            print(tabulate.tabulate(bench[key],
                                    headers=['Change',
                                             f'Before [{hash_1[:8]}] {name_1}',
                                             f'After [{hash_2[:8]}] {name_2}',
                                             'Ratio', 'Benchmark (Parameter)'],
                                    tablefmt="github"))

        return worsened, improved
