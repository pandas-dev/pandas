#!/usr/bin/env python3

# Note: This file has to live next to setup.py or versioneer will not work
import argparse
import os
import sys

import versioneer

sys.path.insert(0, "")


def write_version_info(path) -> None:
    version = None
    git_version = None

    try:
        import _version_meson

        version = _version_meson.__version__
        git_version = _version_meson.__git_version__
    except ImportError:
        version = versioneer.get_version()
        git_version = versioneer.get_versions()["full-revisionid"]
    if os.environ.get("MESON_DIST_ROOT"):
        path = os.path.join(os.environ.get("MESON_DIST_ROOT"), path)
    with open(path, "w", encoding="utf-8") as file:
        file.write(f'__version__="{version}"\n')
        file.write(f'__git_version__="{git_version}"\n')


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        help="Path to write version info to",
        required=False,
    )
    parser.add_argument(
        "--print",
        default=False,
        action="store_true",
        help="Whether to print out the version",
        required=False,
    )
    args = parser.parse_args()

    if args.outfile:
        if not args.outfile.endswith(".py"):
            raise ValueError(
                f"Output file must be a Python file. "
                f"Got: {args.outfile} as filename instead"
            )

        write_version_info(args.outfile)

    if args.print:
        try:
            import _version_meson

            version = _version_meson.__version__
        except ImportError:
            version = versioneer.get_version()
        print(version)


main()
