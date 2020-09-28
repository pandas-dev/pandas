#!/usr/bin/env python3
"""
Tag a pandas release.

This does not push the tag.
"""
import argparse
import subprocess
import sys

from packaging import version


def check_tag(tag):
    assert tag.startswith("v"), "Invalid tag '{}', must start with 'v'".format(tag)
    ver = version.parse(tag.lstrip("v"))
    assert isinstance(ver, version.Version), f"Invalid tag '{tag}'"
    if "rc" in tag:
        assert ".rc" not in tag, "RC tags should be formatted like '.0rcX' "
    return tag


def checkout(tag):
    if tag[-1] == "0" or "rc" in tag:
        # off master
        base = "master"
    else:
        base = ".".join([tag[1:].rsplit(".", 1)[0], "x"])

    subprocess.check_call(["git", "checkout", base])
    # subprocess.check_call(["git", "pull", "--ff-only", "upstream", base])


def commit(tag):
    subprocess.check_call(["git", "clean", "-xdf"])
    print(f"Creating tag {tag}")
    subprocess.check_call(
        [
            "git",
            "commit",
            "--allow-empty",
            "-m",
            "RLS: {}".format(tag[1:]),
        ]
    )
    subprocess.check_call(["git", "tag", "-a", tag, "-m", "Version {}".format(tag[1:])])


def parse_args(args=None):
    parser = argparse.ArgumentParser(__name__, usage=__doc__)
    parser.add_argument("tag", type=check_tag)

    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)
    checkout(args.tag)
    commit(args.tag)


if __name__ == "__main__":
    sys.exit(main())
