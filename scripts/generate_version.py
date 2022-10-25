import argparse
import os

import versioneer


def write_version_info(path):
    if os.environ.get("MESON_DIST_ROOT"):
        # raise ValueError("dist root is", os.environ.get("MESON_DIST_ROOT"))
        path = os.path.join(os.environ.get("MESON_DIST_ROOT"), path)
    with open(path, "w") as file:
        file.write(f'__version__="{versioneer.get_version()}"\n')
        file.write(
            f'__git_version__="{versioneer.get_versions()["full-revisionid"]}"\n'
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o", "--outfile", type=str, help="Path to write version info to"
    )
    args = parser.parse_args()

    if not args.outfile.endswith(".py"):
        raise ValueError(
            f"Output file must be a Python file. "
            f"Got: {args.outfile} as filename instead"
        )

    write_version_info(args.outfile)


main()
