#!/usr/bin/env python3

# Note: This file has to live next to meson.build or versioneer will not work
import argparse
import os
import sys

import versioneer

sys.path.insert(0, "")


def write_version_info(path) -> str:
    version = None
    git_version = None

    try:
        import _version_meson

        version = _version_meson.__version__
        git_version = _version_meson.__git_version__
    except ImportError:
        version = versioneer.get_version()
        git_version = versioneer.get_versions()["full-revisionid"]
    dist_root = os.environ.get("MESON_DIST_ROOT")
    if dist_root:
        path = os.path.join(dist_root, path)
    with open(path, "w", encoding="utf-8") as file:
        file.write(f'__version__="{version}"\n')
        file.write(f'__git_version__="{git_version}"\n')
    return version


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

    version = None
    if args.outfile:
        if not args.outfile.endswith(".py"):
            raise ValueError(
                f"Output file must be a Python file. "
                f"Got: {args.outfile} as filename instead"
            )

        # Reuse the version write_version_info already resolved, rather
        # than re-deriving it below in the `--print` branch. Re-deriving
        # would mean a second (and, in the versioneer fallback path, a
        # third) shell-out to `git describe`/`git rev-parse` when both
        # `--outfile` and `--print` are passed together, for no benefit.
        version = write_version_info(args.outfile)

    if args.print:
        if version is None:
            try:
                import _version_meson

                version = _version_meson.__version__
            except ImportError:
                version = versioneer.get_version()
        print(version)


main()
