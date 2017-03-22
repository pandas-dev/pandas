#!/usr/bin/env python

# script to cancel previous travis builds for the same PR
# originally from
# https://github.com/conda-forge/staged-recipes/pull/2257

try:
    from future_builtins import (
        map,
        filter,
    )
except ImportError:
    pass

import codecs
import contextlib
import json
import os

try:
    from urllib.request import (
        Request,
        urlopen,
    )
except ImportError:
    from urllib2 import (
        Request,
        urlopen,
    )


def check_latest_pr_build(repo, pr, build_num):
    # Not a PR so it is latest.
    if pr is None:
        return True

    headers = {
        "Accept": "application/vnd.travis-ci.2+json",
    }
    url = "https://api.travis-ci.org/repos/{repo}/builds?event_type=pull_request"

    request = Request(url.format(repo=repo), headers=headers)
    with contextlib.closing(urlopen(request)) as response:
        reader = codecs.getreader("utf-8")
        data = json.load(reader(response))

    # Parse the response to get a list of build numbers for this PR.
    builds = data["builds"]
    pr_builds = filter(lambda b: b["pull_request_number"] == pr, builds)
    pr_build_nums = sorted(map(lambda b: int(b["number"]), pr_builds))

    print("build_num: {}".format(build_num))
    print("pr_build_nums: {}".format(','.join([str(n) for n in pr_build_nums])))

    # Check if our build number is the latest (largest)
    # out of all of the builds for this PR.
    if build_num < max(pr_build_nums):
        return False
    else:
        return True


def main():
    repo = os.environ["TRAVIS_REPO_SLUG"]

    pr = os.environ["TRAVIS_PULL_REQUEST"]
    pr = None if pr == "false" else int(pr)
    build_num = int(os.environ["TRAVIS_BUILD_NUMBER"])
    print("checking for fast_finish: {}-{}-{}".format(repo, pr, build_num))

    return int(check_latest_pr_build(repo, pr, build_num) is False)


if __name__ == "__main__":
    import sys
    sys.exit(main())
