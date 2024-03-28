#!/usr/bin/env python3

# This file is vendored from NumPy. See NUMPY_LICENSE.txt

import os
import os.path
import subprocess
import textwrap


def init_version() -> str:
    init = os.path.join(os.path.dirname(__file__), "../../pyproject.toml")
    with open(init, encoding="utf-8") as fid:
        data = fid.readlines()

    version_line = next(line for line in data if line.startswith("version ="))

    result = version_line.strip().split(" = ")[1]
    result = result.replace('"', "").replace("'", "")

    return result


def git_version(version_str: str) -> tuple[str, str]:
    # Append last commit date and hash to dev version information,
    # if available

    ghash = ""
    result = version_str
    try:
        p = subprocess.Popen(
            ["git", "log", "-1", '--format="%H %aI"'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=os.path.dirname(__file__),
        )
    except FileNotFoundError:
        pass
    else:
        out, err = p.communicate()
        if p.returncode == 0:
            ghash, git_date = (
                out.decode("utf-8")
                .strip()
                .replace('"', "")
                .split("T", maxsplit=1)[0]
                .replace("-", "")
                .split()
            )

            # Only attach git tag to development versions
            if "dev" in version_str:
                result = f"{version_str}+git{git_date}.{ghash[:7]}"

    return result, ghash


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--write", help="Save version to this file")
    parser.add_argument(
        "--meson-dist",
        help="Output path is relative to MESON_DIST_ROOT",
        action="store_true",
    )
    args = parser.parse_args()

    version, git_hash = git_version(init_version())

    template = textwrap.dedent(f"""
        version = "{version}"
        __version__ = version
        git_revision = "{git_hash}"
    """)

    if args.write:
        outfile = args.write
        if args.meson_dist:
            outfile = os.path.join(os.environ.get("MESON_DIST_ROOT", ""), outfile)

        # Print human readable output path
        relpath = os.path.relpath(outfile)
        if relpath.startswith("."):
            relpath = outfile

        with open(outfile, "w", encoding="utf-8") as f:
            print(f"Saving version to {relpath}")
            f.write(template)
    else:
        print(version)
