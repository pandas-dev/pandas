from __future__ import annotations

import argparse
import json
import os
import sys

from pathlib import Path
from typing import Any

from setuptools_scm import Configuration
from setuptools_scm._file_finders import find_files
from setuptools_scm._get_version_impl import _get_version
from setuptools_scm.discover import walk_potential_roots


def main(args: list[str] | None = None) -> int:
    opts = _get_cli_opts(args)
    inferred_root: str = opts.root or "."

    pyproject = opts.config or _find_pyproject(inferred_root)

    try:
        config = Configuration.from_file(
            pyproject,
            root=(os.path.abspath(opts.root) if opts.root is not None else None),
        )
    except (LookupError, FileNotFoundError) as ex:
        # no pyproject.toml OR no [tool.setuptools_scm]
        print(
            f"Warning: could not use {os.path.relpath(pyproject)},"
            " using default configuration.\n"
            f" Reason: {ex}.",
            file=sys.stderr,
        )
        config = Configuration(root=inferred_root)
    version: str | None
    if opts.no_version:
        version = "0.0.0+no-version-was-requested.fake-version"
    else:
        version = _get_version(
            config, force_write_version_files=opts.force_write_version_files
        )
    if version is None:
        raise SystemExit("ERROR: no version found for", opts)
    if opts.strip_dev:
        version = version.partition(".dev")[0]

    return command(opts, version, config)


def _get_cli_opts(args: list[str] | None) -> argparse.Namespace:
    prog = "python -m setuptools_scm"
    desc = "Print project version according to SCM metadata"
    parser = argparse.ArgumentParser(prog, description=desc)
    # By default, help for `--help` starts with lower case, so we keep the pattern:
    parser.add_argument(
        "-r",
        "--root",
        default=None,
        help='directory managed by the SCM, default: inferred from config file, or "."',
    )
    parser.add_argument(
        "-c",
        "--config",
        default=None,
        metavar="PATH",
        help="path to 'pyproject.toml' with setuptools-scm config, "
        "default: looked up in the current or parent directories",
    )
    parser.add_argument(
        "--strip-dev",
        action="store_true",
        help="remove the dev/local parts of the version before printing the version",
    )
    parser.add_argument(
        "-N",
        "--no-version",
        action="store_true",
        help="do not include package version in the output",
    )
    output_formats = ["json", "plain", "key-value"]
    parser.add_argument(
        "-f",
        "--format",
        type=str.casefold,
        default="plain",
        help="specify output format",
        choices=output_formats,
    )
    parser.add_argument(
        "-q",
        "--query",
        type=str.casefold,
        nargs="*",
        help="display setuptools-scm settings according to query, "
        "e.g. dist_name, do not supply an argument in order to "
        "print a list of valid queries.",
    )
    parser.add_argument(
        "--force-write-version-files",
        action="store_true",
        help="trigger to write the content of the version files\n"
        "its recommended to use normal/editable installation instead)",
    )
    sub = parser.add_subparsers(title="extra commands", dest="command", metavar="")
    # We avoid `metavar` to prevent printing repetitive information
    desc = "List information about the package, e.g. included files"
    sub.add_parser("ls", help=desc[0].lower() + desc[1:], description=desc)

    # Add create-archival-file subcommand
    archival_desc = "Create .git_archival.txt file for git archive support"
    archival_parser = sub.add_parser(
        "create-archival-file",
        help=archival_desc[0].lower() + archival_desc[1:],
        description=archival_desc,
    )
    archival_group = archival_parser.add_mutually_exclusive_group(required=True)
    archival_group.add_argument(
        "--stable",
        action="store_true",
        help="create stable archival file (recommended, no branch names)",
    )
    archival_group.add_argument(
        "--full",
        action="store_true",
        help="create full archival file with branch information (can cause instability)",
    )
    archival_parser.add_argument(
        "--force", action="store_true", help="overwrite existing .git_archival.txt file"
    )
    return parser.parse_args(args)


# flake8: noqa: C901
def command(opts: argparse.Namespace, version: str, config: Configuration) -> int:
    data: dict[str, Any] = {}

    if opts.command == "ls":
        opts.query = ["files"]

    if opts.command == "create-archival-file":
        return _create_archival_file(opts, config)

    if opts.query == []:
        opts.no_version = True
        sys.stderr.write("Available queries:\n\n")
        opts.query = ["queries"]
        data["queries"] = ["files", *config.__dataclass_fields__]

    if opts.query is None:
        opts.query = []

    if not opts.no_version:
        data["version"] = version

    if "files" in opts.query:
        data["files"] = find_files(config.root)

    for q in opts.query:
        if q in ["files", "queries", "version"]:
            continue

        try:
            if q.startswith("_"):
                raise AttributeError()
            data[q] = getattr(config, q)
        except AttributeError:
            sys.stderr.write(f"Error: unknown query: '{q}'\n")
            return 1

    if opts.format == "json":
        print(json.dumps(data, indent=2))

    if opts.format == "plain":
        _print_plain(data)

    if opts.format == "key-value":
        _print_key_value(data)

    return 0


def _print_plain(data: dict[str, Any]) -> None:
    version = data.pop("version", None)
    if version:
        print(version)
    files = data.pop("files", [])
    for file_ in files:
        print(file_)
    queries = data.pop("queries", [])
    for query in queries:
        print(query)
    if data:
        print("\n".join(data.values()))


def _print_key_value(data: dict[str, Any]) -> None:
    for key, value in data.items():
        if isinstance(value, str):
            print(f"{key} = {value}")
        else:
            str_value = "\n  ".join(value)
            print(f"{key} = {str_value}")


def _find_pyproject(parent: str) -> str:
    for directory in walk_potential_roots(os.path.abspath(parent)):
        pyproject = os.path.join(directory, "pyproject.toml")
        if os.path.isfile(pyproject):
            return pyproject

    return os.path.abspath(
        "pyproject.toml"
    )  # use default name to trigger the default errors


def _create_archival_file(opts: argparse.Namespace, config: Configuration) -> int:
    """Create .git_archival.txt file with appropriate content."""
    archival_path = Path(config.root, ".git_archival.txt")

    # Check if file exists and force flag
    if archival_path.exists() and not opts.force:
        print(
            f"Error: {archival_path} already exists. Use --force to overwrite.",
            file=sys.stderr,
        )
        return 1

    if opts.stable:
        content = _get_stable_archival_content()
        print("Creating stable .git_archival.txt (recommended for releases)")
    elif opts.full:
        content = _get_full_archival_content()
        print("Creating full .git_archival.txt with branch information")
        print("WARNING: This can cause archive checksums to be unstable!")

    try:
        archival_path.write_text(content, encoding="utf-8")
        print(f"Created: {archival_path}")

        gitattributes_path = Path(config.root, ".gitattributes")
        needs_gitattributes = True

        if gitattributes_path.exists():
            # TODO: more nuanced check later
            gitattributes_content = gitattributes_path.read_text("utf-8")
            if (
                ".git_archival.txt" in gitattributes_content
                and "export-subst" in gitattributes_content
            ):
                needs_gitattributes = False

        if needs_gitattributes:
            print("\nNext steps:")
            print("1. Add this line to .gitattributes:")
            print("   .git_archival.txt  export-subst")
            print("2. Commit both files:")
            print("   git add .git_archival.txt .gitattributes")
            print("   git commit -m 'add git archive support'")
        else:
            print("\nNext step:")
            print("Commit the archival file:")
            print("   git add .git_archival.txt")
            print("   git commit -m 'update git archival file'")

        return 0
    except OSError as e:
        print(f"Error: Could not create {archival_path}: {e}", file=sys.stderr)
        return 1


def _get_stable_archival_content() -> str:
    """Generate stable archival file content (no branch names)."""
    return """\
node: $Format:%H$
node-date: $Format:%cI$
describe-name: $Format:%(describe:tags=true,match=*[0-9]*)$
"""


def _get_full_archival_content() -> str:
    """Generate full archival file content with branch information."""
    return """\
# WARNING: Including ref-names can make archive checksums unstable
# after commits are added post-release. Use only if describe-name is insufficient.
node: $Format:%H$
node-date: $Format:%cI$
describe-name: $Format:%(describe:tags=true,match=*[0-9]*)$
ref-names: $Format:%D$
"""
