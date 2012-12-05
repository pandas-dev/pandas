#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

"""Short one-line summary

long summary
"""

import urllib2
import json

import  pandas as pd

def dump_as_gist(data,desc="The Commit"):
    content=dict(version="0.1",timings=data)
    payload=dict(description=desc,
                 public=True,
                 files={'results.txt': dict(content=json.dumps(content))})
    try:
        r=urllib2.urlopen("https://api.github.com/gists", json.dumps(payload))
        if 200 <=r.getcode() <300:
            print("\n\n"+ "-"*80)

            gist = json.loads(r.read())
            file_raw_url = gist['files'].items()[0][1]['raw_url']
            print("[vbench-gist-raw_url] %s" % file_raw_url)
            print("[vbench-html-url] %s" % gist['html_url'])
            print("[vbench-api-url] %s" % gist['url'])

            print("-"*80 +"\n\n")
        else:
            print("api.gothub.com returned status %d" % r.getcode())
    except:
        print("Error occured while dumping to gist")

def main():
    import warnings
    from suite import benchmarks

    warnings.filterwarnings('ignore',category=FutureWarning)

    results=[]
    for b in benchmarks:
        try:
            d=b.run()
            d.update(dict(name=b.name))
            results.append(d)
            msg="{name:<40}: {timing:> 10.4f} [ms]"
            print(msg.format(name=results[-1]['name'],
                             timing=results[-1]['timing']))

        except Exception as e:
            if (type(e) == KeyboardInterrupt or
                'KeyboardInterrupt' in str(d)) :
                raise KeyboardInterrupt()

            msg="{name:<40}: ERROR:\n<-------"
            print(msg.format(name=results[-1]['name']))
            if isinstance(d,dict):
                if d['succeeded']:
                    print("\nException:\n%s\n" % str(e))
                else:
                    for k,v in sorted(d.iteritems()):
                        print("{k}: {v}".format(k=k,v=v))

            print("------->\n")

    dump_as_gist(results,"testing")

def get_vbench_log(build_url):
    r=urllib2.urlopen(build_url)
    if not (200 <= r.getcode() < 300):
        return

    s=json.loads(r.read())
    s=[x for x in s['matrix'] if x['config'].get('env',{}).get('VBENCH')]
            #s=[x for x in s['matrix']]
    if not s:
        return
    id=s[0]['id'] # should be just one for now
    r2=urllib2.urlopen("https://api.travis-ci.org/jobs/%s" % id)
    if (not 200 <= r.getcode() < 300):
        return
    s2=json.loads(r2.read())
    return s2.get('log')

def get_results_raw_url(build):
    import re
    log=get_vbench_log("https://api.travis-ci.org/builds/%s" % build)
    l=[x.strip() for x in log.split("\n") if re.match(".vbench-gist-raw_url",x)]
    if l:
        s=l[0]
        m = re.search("(https://[^\s]+)",s)
        if m:
            return m.group(0)

def get_build_results(build):

    r_url=get_results_raw_url(build)
    if not r_url:
        return
    res=json.loads(urllib2.urlopen(r_url).read())
    timings=res.get("timings")
    if not timings:
        return
    res=[x for x in timings if x.get('succeeded')]
    df = pd.DataFrame(res)
    df = df.set_index("name")
    return df

# builds=[3393087,3393105,3393122,3393125,3393130]
# dfs=[get_build_results(x) for x in builds]
# dfs2=[x[['timing']] for x in dfs]
# for df,b in zip(dfs2,builds):
#     df.columns=[str(b)]
# df = dfs2[0]
# for other in dfs2[1:]:
#     df=df.join(other,how='outer')

if __name__ == "__main__":
	main()
