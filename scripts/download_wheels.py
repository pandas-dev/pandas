#!/usr/bin/env python
"""Fetch wheels from wheels.scipy.org for a pandas version."""
import argparse
import pathlib
import sys
import urllib.parse
import urllib.request

from lxml import html


def parse_args(args=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("version", type=str, help="Pandas version (0.23.0)")
    return parser.parse_args(args)


def fetch(version):
    base = 'http://wheels.scipy.org'
    tree = html.parse(base)
    root = tree.getroot()

    dest = pathlib.Path('dist')
    dest.mkdir(exist_ok=True)

    files = [x for x in root.xpath("//a/text()")
             if x.startswith('pandas-{}'.format(version))
             and not dest.joinpath(x).exists()]

    N = len(files)

    for i, filename in enumerate(files, 1):
        out = str(dest.joinpath(filename))
        link = urllib.request.urljoin(base, filename)
        urllib.request.urlretrieve(link, out)
        print("Downloaded {link} to {out} [{i}/{N}]".format(
            link=link, out=out, i=i, N=N
        ))


def main(args=None):
    args = parse_args(args)
    fetch(args.version)


if __name__ == '__main__':
    sys.exit(main())
