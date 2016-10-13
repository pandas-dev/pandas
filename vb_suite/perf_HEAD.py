#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

"""Run all the vbenches in `suite`, and post the results as a json blob to gist

"""

import urllib2
from contextlib import closing
from urllib2 import urlopen
import json

import pandas as pd

WEB_TIMEOUT = 10


def get_travis_data():
    """figure out what worker we're running on,  and the number of jobs it's running
    """
    import os
    jobid = os.environ.get("TRAVIS_JOB_ID")
    if not jobid:
        return None, None

    with closing(urlopen("https://api.travis-ci.org/workers/")) as resp:
        workers = json.loads(resp.read())

    host = njobs = None
    for item in workers:
        host = item.get("host")
        id = ((item.get("payload") or {}).get("job") or {}).get("id")
        if id and str(id) == str(jobid):
            break
        if host:
            njobs = len(
                [x for x in workers if host in x['host'] and x['payload']])

    return host, njobs


def get_utcdatetime():
    try:
        from datetime import datetime
        return datetime.utcnow().isoformat(" ")
    except:
        pass


def dump_as_gist(data, desc="The Commit", njobs=None):
    host, njobs2 = get_travis_data()[:2]

    if njobs:  # be slightly more reliable
        njobs = max(njobs, njobs2)

    content = dict(version="0.1.1",
                   timings=data,
                   datetime=get_utcdatetime(),   # added in 0.1.1
                   hostname=host,   # added in 0.1.1
                   njobs=njobs    # added in 0.1.1, a measure of load on the travis box
                   )

    payload = dict(description=desc,
                   public=True,
                   files={'results.json': dict(content=json.dumps(content))})
    try:
        with closing(urlopen("https://api.github.com/gists",
                             json.dumps(payload), timeout=WEB_TIMEOUT)) as r:
            if 200 <= r.getcode() < 300:
                print("\n\n" + "-" * 80)

                gist = json.loads(r.read())
                file_raw_url = gist['files'].items()[0][1]['raw_url']
                print("[vbench-gist-raw_url] %s" % file_raw_url)
                print("[vbench-html-url] %s" % gist['html_url'])
                print("[vbench-api-url] %s" % gist['url'])

                print("-" * 80 + "\n\n")
            else:
                print("api.github.com returned status %d" % r.getcode())
    except:
        print("Error occured while dumping to gist")


def main():
    import warnings
    from suite import benchmarks

    exit_code = 0
    warnings.filterwarnings('ignore', category=FutureWarning)

    host, njobs = get_travis_data()[:2]
    results = []
    for b in benchmarks:
        try:
            d = b.run()
            d.update(dict(name=b.name))
            results.append(d)
            msg = "{name:<40}: {timing:> 10.4f} [ms]"
            print(msg.format(name=results[-1]['name'],
                             timing=results[-1]['timing']))

        except Exception as e:
            exit_code = 1
            if (type(e) == KeyboardInterrupt or
                    'KeyboardInterrupt' in str(d)):
                raise KeyboardInterrupt()

            msg = "{name:<40}: ERROR:\n<-------"
            print(msg.format(name=b.name))
            if isinstance(d, dict):
                if d['succeeded']:
                    print("\nException:\n%s\n" % str(e))
                else:
                    for k, v in sorted(d.iteritems()):
                        print("{k}: {v}".format(k=k, v=v))

            print("------->\n")

    dump_as_gist(results, "testing", njobs=njobs)

    return exit_code


if __name__ == "__main__":
    import sys
    sys.exit(main())

#####################################################
# functions for retrieving and processing the results


def get_vbench_log(build_url):
    with closing(urllib2.urlopen(build_url)) as r:
        if not (200 <= r.getcode() < 300):
            return

        s = json.loads(r.read())
        s = [x for x in s['matrix'] if "VBENCH" in ((x.get('config', {})
                                                    or {}).get('env', {}) or {})]
                # s=[x for x in s['matrix']]
        if not s:
            return
        id = s[0]['id']  # should be just one for now
        with closing(urllib2.urlopen("https://api.travis-ci.org/jobs/%s" % id)) as r2:
            if not 200 <= r.getcode() < 300:
                return
            s2 = json.loads(r2.read())
            return s2.get('log')


def get_results_raw_url(build):
    "Taks a Travis a build number, retrieves the build log and extracts the gist url"
    import re
    log = get_vbench_log("https://api.travis-ci.org/builds/%s" % build)
    if not log:
        return
    l = [x.strip(
    ) for x in log.split("\n") if re.match(".vbench-gist-raw_url", x)]
    if l:
        s = l[0]
        m = re.search("(https://[^\s]+)", s)
        if m:
            return m.group(0)


def convert_json_to_df(results_url):
    """retrieve json results file from url and return df

    df contains timings for all successful vbenchmarks
    """

    with closing(urlopen(results_url)) as resp:
        res = json.loads(resp.read())
    timings = res.get("timings")
    if not timings:
        return
    res = [x for x in timings if x.get('succeeded')]
    df = pd.DataFrame(res)
    df = df.set_index("name")
    return df


def get_build_results(build):
    "Returns a df with the results of the VBENCH job associated with the travis build"
    r_url = get_results_raw_url(build)
    if not r_url:
        return

    return convert_json_to_df(r_url)


def get_all_results(repo_id=53976):  # travis pandas-dev/pandas id
    """Fetches the VBENCH results for all travis builds, and returns a list of result df

   unsuccesful individual vbenches are dropped.
    """
    from collections import OrderedDict

    def get_results_from_builds(builds):
        dfs = OrderedDict()
        for build in builds:
            build_id = build['id']
            build_number = build['number']
            print(build_number)
            res = get_build_results(build_id)
            if res is not None:
                dfs[build_number] = res
        return dfs

    base_url = 'https://api.travis-ci.org/builds?url=%2Fbuilds&repository_id={repo_id}'
    url = base_url.format(repo_id=repo_id)
    url_after = url + '&after_number={after}'
    dfs = OrderedDict()

    while True:
        with closing(urlopen(url)) as r:
            if not (200 <= r.getcode() < 300):
                break
            builds = json.loads(r.read())
        res = get_results_from_builds(builds)
        if not res:
            break
        last_build_number = min(res.keys())
        dfs.update(res)
        url = url_after.format(after=last_build_number)

    return dfs


def get_all_results_joined(repo_id=53976):
    def mk_unique(df):
        for dupe in df.index.get_duplicates():
            df = df.ix[df.index != dupe]
        return df
    dfs = get_all_results(repo_id)
    for k in dfs:
        dfs[k] = mk_unique(dfs[k])
    ss = [pd.Series(v.timing, name=k) for k, v in dfs.iteritems()]
    results = pd.concat(reversed(ss), 1)
    return results
