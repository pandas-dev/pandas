import argparse
import os
import sys
import shutil
import subprocess
from typing import Iterable, Iterator, List, Optional, Text, Tuple

from .color import color_unified_diff_line
from .diff import run_external_diff, u_diff
from .utils import file_exists, get_tables_argument_list


def pipe_output(output: str) -> None:
    """Pipes output to a pager if stdout is a TTY and a pager is available."""

    if not output:
        return

    if not sys.stdout.isatty():
        sys.stdout.write(output)
        return

    pager = os.getenv("PAGER") or shutil.which("less")

    if not pager:
        sys.stdout.write(output)
        return

    pager_cmd = [pager]
    if "less" in os.path.basename(pager):
        pager_cmd.append("-R")

    proc = subprocess.Popen(pager_cmd, stdin=subprocess.PIPE, text=True)
    try:
        proc.stdin.write(output)
        proc.stdin.close()
        proc.wait()
    except (BrokenPipeError, KeyboardInterrupt):
        # Pager process was terminated before all output was written.
        # This is not an error. The main exception handler will deal with it.
        if proc.stdin:
            proc.stdin.close()
        # The process might still be running, but we have closed our side of the
        # pipe. The Popen destructor will send a SIGKILL to the child.
    except Exception:
        if proc.stdin:
            proc.stdin.close()
        raise


def _is_gnu_diff(diff_tool: str) -> bool:
    """Returns True if the provided diff executable is GNU diff."""
    try:
        proc = subprocess.run(
            [diff_tool, "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except OSError:
        return False

    version_output = (proc.stdout or "") + (proc.stderr or "")
    return "GNU diffutils" in version_output


def _iter_filtered_table_tags(
    tags: Iterable[str],
    include_tables: Optional[List[str]] = None,
    exclude_tables: Optional[List[str]] = None,
) -> Iterator[str]:
    for tag in tags:
        if exclude_tables and tag in exclude_tables:
            continue
        if include_tables and tag not in include_tables:
            continue
        yield tag


def summarize(
    file1: str,
    file2: str,
    include_tables: Optional[List[str]] = None,
    exclude_tables: Optional[List[str]] = None,
    font_number_1: int = -1,
    font_number_2: int = -1,
) -> Tuple[bool, str]:
    from fontTools.ttLib import TTFont

    with (
        TTFont(file1, lazy=True, fontNumber=font_number_1) as font1,
        TTFont(file2, lazy=True, fontNumber=font_number_2) as font2,
    ):
        tags1 = {str(tag) for tag in font1.reader.keys()}
        tags2 = {str(tag) for tag in font2.reader.keys()}

        all_tags = sorted(
            set(
                _iter_filtered_table_tags(
                    tags1 | tags2,
                    include_tables=include_tables,
                    exclude_tables=exclude_tables,
                )
            )
        )

        only1 = [tag for tag in all_tags if tag in tags1 and tag not in tags2]
        only2 = [tag for tag in all_tags if tag in tags2 and tag not in tags1]
        both = [tag for tag in all_tags if tag in tags1 and tag in tags2]

        identical = True
        lines: List[str] = []

        lines.append(f"Binary table summary:\n")
        lines.append(f"  file1: {file1}\n")
        lines.append(f"  file2: {file2}\n")

        if only1:
            identical = False
            lines.append(f"\nTables only in file1 ({len(only1)}):\n")
            for tag in only1:
                lines.append(f"- {tag} ({len(font1.reader[tag])} bytes)\n")
        if only2:
            identical = False
            lines.append(f"\nTables only in file2 ({len(only2)}):\n")
            for tag in only2:
                lines.append(f"+ {tag} ({len(font2.reader[tag])} bytes)\n")

        lines.append(f"\nTables in both ({len(both)}):\n")
        for tag in both:
            data1 = font1.reader[tag]
            data2 = font2.reader[tag]
            if data1 == data2:
                lines.append(f"  {tag}: SAME ({len(data1)} bytes)\n")
            else:
                identical = False
                lines.append(f"* {tag}: DIFF ({len(data1)} vs {len(data2)} bytes)\n")

        if identical:
            lines.append("\nResult: SAME\n")
        else:
            lines.append("\nResult: DIFFERENT\n")

        return identical, "".join(lines)


def get_binary_exclude_tables(
    file1: str,
    file2: str,
    include_tables: Optional[List[str]] = None,
    exclude_tables: Optional[List[str]] = None,
    font_number_1: int = -1,
    font_number_2: int = -1,
) -> Tuple[bool, str]:
    from fontTools.ttLib import TTFont

    with (
        TTFont(file1, lazy=True, fontNumber=font_number_1) as font1,
        TTFont(file2, lazy=True, fontNumber=font_number_2) as font2,
    ):
        tags1 = {str(tag) for tag in font1.reader.keys()}
        tags2 = {str(tag) for tag in font2.reader.keys()}

        all_tags = sorted(
            set(
                _iter_filtered_table_tags(
                    tags1 | tags2,
                    include_tables=include_tables,
                    exclude_tables=exclude_tables,
                )
            )
        )

        both = [tag for tag in all_tags if tag in tags1 and tag in tags2]
        out = set()

        for tag in both:
            data1 = font1.reader[tag]
            data2 = font2.reader[tag]
            if data1 == data2:
                out.add(tag)

        return out


def main():
    """Compare two fonts for differences"""
    # try/except block rationale:
    # handles "premature" socket closure exception that is
    # raised by Python when stdout is piped to tools like
    # the `head` executable and socket is closed early
    # see: https://docs.python.org/3/library/signal.html#note-on-sigpipe
    ret = 0
    try:
        ret = run(sys.argv[1:])
    except KeyboardInterrupt:
        pass
    except BrokenPipeError:
        # Python flushes standard streams on exit; redirect remaining output
        # to devnull to avoid another BrokenPipeError at shutdown
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
    return ret


def run(argv: List[Text]):
    # ------------------------------------------
    # argparse command line argument definitions
    # ------------------------------------------
    parser = argparse.ArgumentParser(
        description="An OpenType table diff tool for fonts."
    )
    parser.add_argument(
        "-l",
        "--summary",
        action="store_true",
        help="Report table presence and binary equality only",
    )
    parser.add_argument(
        "-U",
        "--lines",
        type=int,
        default=3,
        help="Number of context lines for unified diff (default: 3)",
    )
    parser.add_argument(
        "-t",
        "--include",
        type=str,
        nargs="+",
        default=None,
        help="Font tables to include. Multiple options are allowed.",
    )
    parser.add_argument(
        "-x",
        "--exclude",
        type=str,
        nargs="+",
        default=None,
        help="Font tables to exclude. Multiple options are allowed.",
    )
    parser.add_argument(
        "--diff", type=str, help="Run external diff tool command (default: diff)"
    )
    parser.add_argument(
        "--diff-arg",
        type=str,
        default=None,
        help="External diff tool arguments (default: -u)",
    )
    parser.add_argument(
        "--color",
        choices=["auto", "never", "always"],
        default="auto",
        help="Whether to colorize output (default: auto)",
    )
    parser.add_argument(
        "--y1",
        type=int,
        default=-1,
        metavar="NUMBER",
        help="Select font number for TrueType Collection (.ttc/.otc) FILE1, starting from 0",
    )
    parser.add_argument(
        "--y2",
        type=int,
        default=-1,
        metavar="NUMBER",
        help="Select font number for TrueType Collection (.ttc/.otc) FILE2, starting from 0",
    )
    parser.add_argument(
        "-a",
        "--always",
        action="store_true",
        help="Compare tables even if binary identical",
    )
    parser.add_argument(
        "-b",
        "--binary",
        action="store_true",
        help="Compare tables only if binaries differ (default)",
    )
    parser.add_argument(
        "-q", "--quiet", action="store_true", help="Suppress all output"
    )
    parser.add_argument("FILE1", help="Font file path 1")
    parser.add_argument("FILE2", help="Font file path 2")

    args: argparse.Namespace = parser.parse_args(argv)

    # /////////////////////////////////////////////////////////
    #
    #  Validations
    #
    # /////////////////////////////////////////////////////////

    # ----------------------------------
    #  Incompatible argument validations
    # ----------------------------------

    if args.always and args.binary:
        if not args.quiet:
            sys.stderr.write(
                f"[*] Error: --always and --binary are mutually exclusive options. "
                f"Please use ONLY one of these options in your command.{os.linesep}"
            )
        return 2
    if not args.always:
        args.binary = True

    # -------------------------------
    #  File path argument validations
    # -------------------------------

    if not file_exists(args.FILE1):
        if not args.quiet:
            sys.stderr.write(
                f"[*] ERROR: The file path '{args.FILE1}' can not be found.{os.linesep}"
            )
        return 2
    if not file_exists(args.FILE2):
        if not args.quiet:
            sys.stderr.write(
                f"[*] ERROR: The file path '{args.FILE2}' can not be found.{os.linesep}"
            )
        return 2

    # /////////////////////////////////////////////////////////
    #
    #  Command line logic
    #
    # /////////////////////////////////////////////////////////

    # parse explicitly included or excluded tables in
    # the command line arguments
    # set as a Python list if it was defined on the command line
    # or as None if it was not set on the command line
    include_list: Optional[List[Text]] = get_tables_argument_list(args.include)
    exclude_list: Optional[List[Text]] = get_tables_argument_list(args.exclude)

    if args.summary:
        try:
            identical, output = summarize(
                args.FILE1,
                args.FILE2,
                include_tables=include_list,
                exclude_tables=exclude_list,
                font_number_1=args.y1,
                font_number_2=args.y2,
            )
            if not args.quiet:
                sys.stdout.write(output)
            return 0 if identical else 1
        except Exception as e:
            if not args.quiet:
                sys.stderr.write(f"[*] ERROR: {e}{os.linesep}")
            return 2

    if args.binary:
        excluded_binary_tables = get_binary_exclude_tables(
            args.FILE1,
            args.FILE2,
            include_tables=include_list,
            exclude_tables=exclude_list,
            font_number_1=args.y1,
            font_number_2=args.y2,
        )
        if include_list is not None:
            include_list = [
                tag for tag in include_list if tag not in excluded_binary_tables
            ]
        else:
            if exclude_list is None:
                exclude_list = []
            exclude_list.extend(sorted(excluded_binary_tables))

    diff_tool = args.diff
    color_output = args.color == "always" or (
        args.color == "auto" and sys.stdout.isatty
    )

    if diff_tool is None:
        diff_tool = shutil.which("diff")
    elif diff_tool:
        diff_tool = shutil.which(diff_tool)
        if diff_tool is None:
            if not args.quiet:
                sys.stderr.write(
                    f"[*] ERROR: The external diff tool executable "
                    f"'{args.diff}' was not found.{os.linesep}"
                )
            return 2

    try:
        if diff_tool:
            diff_arg = args.diff_arg
            if diff_arg is None:
                if args.lines == 3:
                    diff_arg = ["-u"]
                else:
                    diff_arg = ["-u{}".format(args.lines)]
                if _is_gnu_diff(diff_tool):
                    diff_arg.append(r"-F^\s\s<")
            else:
                diff_arg = diff_arg.split()

            output = run_external_diff(
                diff_tool,
                diff_arg,
                args.FILE1,
                args.FILE2,
                include_tables=include_list,
                exclude_tables=exclude_list,
                font_number_a=args.y1,
                font_number_b=args.y2,
                use_multiprocess=True,
            )
        else:
            output = u_diff(
                args.FILE1,
                args.FILE2,
                context_lines=args.lines,
                include_tables=include_list,
                exclude_tables=exclude_list,
                font_number_a=args.y1,
                font_number_b=args.y2,
                use_multiprocess=True,
            )

        if color_output:
            output = [color_unified_diff_line(line) for line in output]

        output = "".join(output)
        if not args.quiet:
            pipe_output(output)
        return 1 if output else 0

    except Exception as e:
        if not args.quiet:
            sys.stderr.write(f"[*] ERROR: {e}{os.linesep}")
        return 2
