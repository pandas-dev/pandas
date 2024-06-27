# Licensed under a 3-clause BSD style license FOR Mamba - see LICENSE in mamba-org/mamba
# Also covered by the 3-clause BSD style license FOR boa - see LICENSE in mamba-org/boa
# Copyright 2019 QuantStack and the Mamba contributors.
# Very lightly edited / simplified for use within asv
import os
import urllib.parse
import collections

import libmambapy
from conda.base.constants import ChannelPriority
from conda.base.context import context
from conda.core.index import check_allowlist
from conda.gateways.connection.session import CondaHttpAuth


def get_index(
    channel_urls=(),
    prepend=True,
    platform=None,
    use_local=False,
    use_cache=False,
    unknown=None,
    prefix=None,
    repodata_fn="repodata.json",
):
    if isinstance(platform, str):
        platform = [platform, "noarch"]

    all_channels = []
    if use_local:
        all_channels.append("local")
    all_channels.extend(channel_urls)
    if prepend:
        all_channels.extend(context.channels)
    check_allowlist(all_channels)

    # Remove duplicates but retain order
    all_channels = list(collections.OrderedDict.fromkeys(all_channels))

    dlist = libmambapy.DownloadTargetList()

    index = []

    def fixup_channel_spec(spec):
        at_count = spec.count("@")
        if at_count > 1:
            first_at = spec.find("@")
            spec = (
                spec[:first_at]
                + urllib.parse.quote(spec[first_at])
                + spec[first_at + 1 :]
            )
        if platform:
            spec = spec + "[" + ",".join(platform) + "]"
        return spec

    all_channels = list(map(fixup_channel_spec, all_channels))
    pkgs_dirs = libmambapy.MultiPackageCache(context.pkgs_dirs)
    libmambapy.create_cache_dir(str(pkgs_dirs.first_writable_path))

    for channel in libmambapy.get_channels(all_channels):
        for channel_platform, url in channel.platform_urls(with_credentials=True):
            full_url = CondaHttpAuth.add_binstar_token(url)

            sd = libmambapy.SubdirData(
                channel, channel_platform, full_url, pkgs_dirs, repodata_fn
            )

            needs_finalising = sd.download_and_check_targets(dlist)
            index.append(
                (
                    sd,
                    {
                        "platform": channel_platform,
                        "url": url,
                        "channel": channel,
                        "needs_finalising": needs_finalising,
                    },
                )
            )

    for sd, info in index:
        if info["needs_finalising"]:
            sd.finalize_checks()
        dlist.add(sd)

    is_downloaded = dlist.download(libmambapy.MAMBA_DOWNLOAD_FAILFAST)

    if not is_downloaded:
        raise RuntimeError("Error downloading repodata.")

    return index


def load_channels(
    pool,
    channels,
    repos,
    has_priority=None,
    prepend=True,
    platform=None,
    use_local=False,
    use_cache=True,
    repodata_fn="repodata.json",
):
    index = get_index(
        channel_urls=channels,
        prepend=prepend,
        platform=platform,
        use_local=use_local,
        repodata_fn=repodata_fn,
        use_cache=use_cache,
    )

    if has_priority is None:
        has_priority = context.channel_priority in [
            ChannelPriority.STRICT,
            ChannelPriority.FLEXIBLE,
        ]

    subprio_index = len(index)
    if has_priority:
        # first, count unique channels
        n_channels = len(set([entry["channel"].canonical_name for _, entry in index]))
        current_channel = index[0][1]["channel"].canonical_name
        channel_prio = n_channels

    for subdir, entry in index:
        # add priority here
        if has_priority:
            if entry["channel"].canonical_name != current_channel:
                channel_prio -= 1
                current_channel = entry["channel"].canonical_name
            priority = channel_prio
        else:
            priority = 0
        if has_priority:
            subpriority = 0
        else:
            subpriority = subprio_index
            subprio_index -= 1

        if not subdir.loaded() and entry["platform"] != "noarch":
            # ignore non-loaded subdir if channel is != noarch
            continue

        if context.verbosity != 0 and not context.json:
            print(
                "Channel: {}, platform: {}, prio: {} : {}".format(
                    entry["channel"], entry["platform"], priority, subpriority
                )
            )
            print("Cache path: ", subdir.cache_path())

        repo = subdir.create_repo(pool)
        repo.set_priority(priority, subpriority)
        repos.append(repo)

    return index


class MambaSolver:
    def __init__(self, channels, platform, context, output_folder=None):
        self.channels = channels
        self.platform = platform
        self.context = context
        self.output_folder = output_folder or "local"
        self.pool = libmambapy.Pool()
        self.repos = []
        self.index = load_channels(
            self.pool, self.channels, self.repos, platform=platform
        )

        self.local_index = []
        self.local_repos = {}
        # load local repo, too
        self.replace_channels()

    def replace_installed(self, prefix):
        prefix_data = libmambapy.PrefixData(prefix)
        vp = libmambapy.get_virtual_packages()
        prefix_data.add_packages(vp)
        repo = libmambapy.Repo(self.pool, prefix_data)
        repo.set_installed()

    def replace_channels(self):
        self.local_index = get_index(
            (self.output_folder,), platform=self.platform, prepend=False
        )

        for _, v in self.local_repos.items():
            v.clear(True)

        start_prio = len(self.channels) + len(self.index)
        for subdir, channel in self.local_index:
            if not subdir.loaded():
                continue

            # support new mamba
            if isinstance(channel, dict):
                channelstr = channel["url"]
                channelurl = channel["url"]
            else:
                channelstr = str(channel)
                channelurl = channel.url(with_credentials=True)

            cp = subdir.cache_path()
            if cp.endswith(".solv"):
                os.remove(subdir.cache_path())
                cp = cp.replace(".solv", ".json")

            self.local_repos[channelstr] = libmambapy.Repo(
                self.pool, channelstr, cp, channelurl
            )

            self.local_repos[channelstr].set_priority(start_prio, 0)
            start_prio -= 1

    def solve(self, specs, pkg_cache_path=None):
        """Solve given a set of specs.
        Parameters
        ----------
        specs : list of str
            A list of package specs. You can use `conda.models.match_spec.MatchSpec`
            to get them to the right form by calling
            `MatchSpec(mypec).conda_build_form()`
        Returns
        -------
        transaction : libmambapy.Transaction
            The mamba transaction.
        Raises
        ------
        RuntimeError :
            If the solver did not find a solution.
        """
        solver_options = [(libmambapy.SOLVER_FLAG_ALLOW_DOWNGRADE, 1)]
        api_solver = libmambapy.Solver(self.pool, solver_options)
        _specs = specs

        api_solver.add_jobs(_specs, libmambapy.SOLVER_INSTALL)
        success = api_solver.try_solve()

        if not success:
            error_string = "Mamba failed to solve:\n"
            for s in _specs:
                error_string += f" - {s}\n"
            error_string += "\nwith channels:\n"
            for c in self.channels:
                error_string += f" - {c}\n"
            error_string += api_solver.explain_problems()
            print(error_string)
            raise RuntimeError("Solver could not find solution." + error_string)

        if pkg_cache_path is None:
            # use values from conda
            pkg_cache_path = self.context.pkgs_dirs

        package_cache = libmambapy.MultiPackageCache(pkg_cache_path)
        return libmambapy.Transaction(api_solver, package_cache)
