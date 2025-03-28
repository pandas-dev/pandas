# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import re
import itertools
import datetime
import urllib.parse

from ..results import iter_results
from ..console import log
from ..publishing import OutputPublisher
from ..step_detect import detect_regressions
from .. import util, feed


class Regressions(OutputPublisher):
    name = "regressions"
    button_label = "Show regressions"
    description = "Display information about recent regressions"
    order = 3

    @classmethod
    def publish(cls, conf, repo, benchmarks, graphs, revisions):
        # Analyze the data in the graphs --- it's been cleaned up and
        # it's easier to work with than the results directly

        regressions = []
        revision_to_hash = dict((r, h) for h, r in revisions.items())

        data_filter = _GraphDataFilter(conf, repo, revisions)

        all_params = graphs.get_params()

        for j, (file_name, graph) in enumerate(graphs):
            if 'summary' in graph.params:
                continue

            benchmark_name = os.path.basename(file_name)
            benchmark = benchmarks.get(benchmark_name)
            if not benchmark:
                continue

            log.dot()

            for graph_data in data_filter.get_graph_data(graph, benchmark):
                cls._process_regression(regressions, revision_to_hash, repo, all_params,
                                        graph_data, graph)

        cls._save(conf, {'regressions': regressions})
        cls._save_feed(conf, benchmarks, regressions, revisions, revision_to_hash)

    @classmethod
    def _process_regression(cls, regressions, revision_to_hash, repo,
                            all_params, graph_data, graph):
        j, entry_name, steps, threshold = graph_data

        last_v, best_v, jumps = detect_regressions(steps, threshold)

        if last_v is None:
            return

        # Select unique graph params
        graph_params = {}
        for name, value in graph.params.items():
            if len(all_params[name]) > 1:
                graph_params[name] = value

        graph_path = graph.path + '.json'

        # Check which ranges are a single commit
        for k, jump in enumerate(jumps):
            commit_a = revision_to_hash[jump[0]]
            commit_b = revision_to_hash[jump[1]]
            spec = repo.get_range_spec(commit_a, commit_b)
            commits = repo.get_hashes_from_range(spec)
            if len(commits) == 1:
                jumps[k] = (None, jump[1], jump[2], jump[3])

        # Produce output
        regression = [entry_name, graph_path, graph_params, j, last_v, best_v, jumps]
        regressions.append(regression)

    @classmethod
    def _save(cls, conf, data):
        fn = os.path.join(conf.html_dir, 'regressions.json')
        util.write_json(fn, data, compact=True)

    @classmethod
    def _save_feed(cls, conf, benchmarks, data, revisions, revision_to_hash):
        """
        Save the results as an Atom feed
        """

        filename = os.path.join(conf.html_dir, 'regressions.xml')

        # Determine publication date as the date when the benchmark
        # was run --- if it is missing, use the date of the commit
        run_timestamps = {}
        revision_timestamps = {}
        for results in iter_results(conf.results_dir):
            if results.commit_hash not in revisions:
                # revisions could be filtered when specifying a range
                # in 'asv publish'
                continue
            revision = revisions[results.commit_hash]
            revision_timestamps[revision] = results.date

            # Time when the benchmark was run
            for benchmark_name, timestamp in results.started_at.items():
                if timestamp is None:
                    continue
                key = (benchmark_name, revision)
                run_timestamps[key] = timestamp

            # Fallback to commit date
            for benchmark_name in results.get_result_keys(benchmarks):
                key = (benchmark_name, revision)
                run_timestamps.setdefault(key, results.date)

        # Generate feed entries
        entries = []

        for name, graph_path, graph_params, idx, last_value, best_value, jumps in data:
            if '(' in name:
                benchmark_name = name[:name.index('(')]
            else:
                benchmark_name = name

            benchmark = benchmarks[benchmark_name]

            if idx is not None:
                graph_params = dict(graph_params)

                # Add URL parameters
                param_values, = itertools.islice(itertools.product(*benchmark['params']),
                                                 idx, idx + 1)
                for k, v in zip(benchmark['param_names'], param_values):
                    graph_params['p-' + k] = v

            for rev1, rev2, value1, value2 in jumps:
                timestamps = (run_timestamps[benchmark_name, t]
                              for t in (rev1, rev2) if t is not None)
                last_timestamp = max(timestamps)

                updated = datetime.datetime.fromtimestamp(last_timestamp / 1000)

                params = dict(graph_params)

                if rev1 is None:
                    params['commits'] = f'{revision_to_hash[rev2]}'
                else:
                    params['commits'] = '{0}-{1}'.format(revision_to_hash[rev1],
                                                         revision_to_hash[rev2])

                link = f'index.html#{benchmark_name}?{urllib.parse.urlencode(params)}'

                try:
                    best_percentage = "{0:.2f}%".format(100 *
                                                        (last_value - best_value) / best_value)
                except ZeroDivisionError:
                    best_percentage = f"{last_value - best_value:.2g} units"

                try:
                    percentage = f"{100 * (value2 - value1) / value1:.2f}%"
                except ZeroDivisionError:
                    percentage = f"{value2 - value1:.2g} units"

                jump_date = datetime.datetime.fromtimestamp(revision_timestamps[rev2] / 1000)
                jump_date_str = jump_date.strftime('%Y-%m-%d %H:%M:%S')

                if rev1 is not None:
                    commit_a = revision_to_hash[rev1]
                    commit_b = revision_to_hash[rev2]
                    if 'github.com' in conf.show_commit_url:
                        commit_url = (conf.show_commit_url + '../compare/' +
                                      commit_a + "..." + commit_b)
                    else:
                        commit_url = conf.show_commit_url + commit_a
                    commit_ref = 'in commits <a href="{0}">{1}...{2}</a>'.format(commit_url,
                                                                                 commit_a[:8],
                                                                                 commit_b[:8])
                else:
                    commit_a = revision_to_hash[rev2]
                    commit_url = conf.show_commit_url + commit_a
                    commit_ref = f'in commit <a href="{commit_url}">{commit_a[:8]}</a>'

                unit = benchmark.get('unit', '')
                best_value_str = util.human_value(best_value, unit)
                last_value_str = util.human_value(last_value, unit)
                value1_str = util.human_value(value1, unit)
                value2_str = util.human_value(value2, unit)

                title = "{percentage} {name}".format(**locals())
                summary = """
                <a href="{link}">{percentage} regression</a> on {jump_date_str} {commit_ref}.<br>
                New value: {value2_str}, old value: {value1_str}.<br>
                Latest value: {last_value_str} ({best_percentage} worse
                than best value {best_value_str}).
                """.format(**locals()).strip()

                # Information that uniquely identifies a regression
                # --- if the values and the texts change on later
                # runs, feed readers should is identify the regression
                # as the same one, as long as the benchmark name and
                # commits match.
                id_context = [name, revision_to_hash.get(rev1, ""), revision_to_hash.get(rev2, "")]
                id_date = util.js_timestamp_to_datetime(revision_timestamps[rev2])

                entries.append(feed.FeedEntry(title, updated, link, summary, id_context, id_date))

        entries.sort(key=lambda x: x.updated, reverse=True)

        feed.write_atom(filename, entries,
                        title=f'{conf.project} performance regressions',
                        author='Airspeed Velocity',
                        address=f'{conf.project}.asv')


class _GraphDataFilter:
    """
    Obtain data sets from graphs, following configuration settings.
    """

    def __init__(self, conf, repo, revisions):
        self.conf = conf
        self.repo = repo
        self.revisions = revisions
        self._start_revisions = {}

    def get_graph_data(self, graph, benchmark):
        """
        Iterator over graph data sets

        Yields
        ------
        param_idx
            Flat index to parameter permutations for parameterized benchmarks.
            None if benchmark is not parameterized.
        entry_name
            Name for the data set. If benchmark is non-parameterized, this is the
            benchmark name.
        steps
            Steps to consider in regression detection.
        threshold
            User-specified threshold for regression detection.

        """
        if benchmark.get('params'):
            param_iter = enumerate(zip(itertools.product(*benchmark['params']),
                                       graph.get_steps()))
        else:
            param_iter = [(None, (None, graph.get_steps()))]

        for j, (param, steps) in param_iter:
            if param is None:
                entry_name = benchmark['name']
            else:
                entry_name = benchmark['name'] + f"({', '.join(param)})"

            start_revision = self._get_start_revision(graph, benchmark, entry_name)
            threshold = self._get_threshold(graph, benchmark, entry_name)

            if start_revision is None:
                # Skip detection
                continue

            steps = [step for step in steps if step[1] >= start_revision]

            yield j, entry_name, steps, threshold

    def _get_start_revision(self, graph, benchmark, entry_name):
        """
        Compute the first revision allowed by asv.conf.json.

        Revisions correspond to linearized commit history and the
        regression detection runs on this order --- the starting commit
        thus corresponds to a specific starting revision.
        """
        start_revision = min(self.revisions.values())

        if graph.params.get('branch'):
            branch_suffix = '@' + graph.params.get('branch')
        else:
            branch_suffix = ''

        for regex, start_commit in self.conf.regressions_first_commits.items():
            if re.match(regex, entry_name + branch_suffix):
                if start_commit is None:
                    # Disable regression detection completely
                    return None

                if self.conf.branches == [None]:
                    key = (start_commit, None)
                else:
                    key = (start_commit, graph.params.get('branch'))

                if key not in self._start_revisions:
                    spec = self.repo.get_new_range_spec(*key)
                    start_hash = self.repo.get_hash_from_name(start_commit)

                    for commit in [start_hash] + self.repo.get_hashes_from_range(spec):
                        rev = self.revisions.get(commit)
                        if rev is not None:
                            self._start_revisions[key] = rev
                            break
                    else:
                        # Commit not found in the branch --- warn and ignore.
                        log.warning(("Commit {0} specified in `regressions_first_commits` "
                                     "not found in branch").format(start_commit))
                        self._start_revisions[key] = -1

                start_revision = max(start_revision, self._start_revisions[key] + 1)

        return start_revision

    def _get_threshold(self, graph, benchmark, entry_name):
        """
        Compute the regression threshold in asv.conf.json.
        """
        if graph.params.get('branch'):
            branch_suffix = '@' + graph.params.get('branch')
        else:
            branch_suffix = ''

        max_threshold = None

        for regex, threshold in self.conf.regressions_thresholds.items():
            if re.match(regex, entry_name + branch_suffix):
                try:
                    threshold = float(threshold)
                except ValueError:
                    raise util.UserError("Non-float threshold in asv.conf.json: {!r}"
                                         .format(threshold))

                if max_threshold is None:
                    max_threshold = threshold
                else:
                    max_threshold = max(threshold, max_threshold)

        if max_threshold is None:
            max_threshold = 0.05

        return max_threshold
