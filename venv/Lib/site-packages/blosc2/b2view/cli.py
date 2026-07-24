"""Command line entry point for b2view."""

from __future__ import annotations

import argparse
import os
import sys

# Demo server roots: --download takes a path relative to the @public root (e.g.
# "large/chicago-taxi-flat.b2z").  The info endpoint returns the bundle's
# metadata (notably "cbytes", the stored size), giving a determinate progress bar.
DEFAULT_DOWNLOAD_PATH = "large/chicago-taxi-flat.b2z"
DOWNLOAD_BASE_URL = "https://cat2.cloud/demo/api/download/@public/"
INFO_BASE_URL = "https://cat2.cloud/demo/api/info/@public/"


def resolve_source(
    urlpath: str | None, download: str | None, *, exists=os.path.exists
) -> tuple[str, str | None, str | None]:
    """Resolve the CLI args to ``(urlpath, download_url, info_url)``.

    Without ``--download`` the positional *urlpath* is used as-is.  With
    ``--download PATH`` (a path relative to the demo's ``@public`` root), the
    bundle is saved in the cwd under its basename; if that file is already there
    the download is skipped (both URLs None), otherwise it is fetched from
    :data:`DOWNLOAD_BASE_URL` and its size read from :data:`INFO_BASE_URL`.  The
    two arguments are mutually exclusive.
    """
    if download is not None:
        if urlpath is not None:
            raise ValueError("--download cannot be combined with a positional path")
        dest = os.path.basename(download)  # local file lands in the cwd
        if exists(dest):
            return dest, None, None
        return dest, DOWNLOAD_BASE_URL + download, INFO_BASE_URL + download
    if urlpath is None:
        raise ValueError("provide a path to a .b2d/.b2z bundle, or use --download")
    return urlpath, None, None


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Browse a Blosc2 TreeStore bundle in the terminal.")
    parser.add_argument("urlpath", nargs="?", default=None, help="Path to a .b2d directory or .b2z file")
    parser.add_argument("path", nargs="?", default="/", help="Optional starting path inside the bundle")
    parser.add_argument(
        "--download",
        nargs="?",
        const=DEFAULT_DOWNLOAD_PATH,
        default=None,
        metavar="PATH",
        help=(
            f"Open the bundle at PATH (relative to the demo server's @public root; "
            f"default {DEFAULT_DOWNLOAD_PATH!r}) from the cwd, downloading it first if not present"
        ),
    )
    parser.add_argument("--preview-rows", type=int, default=20, help="Maximum preview rows")
    parser.add_argument("--preview-cols", type=int, default=10, help="Maximum preview columns")
    parser.add_argument(
        "--panel",
        choices=["tree", "meta", "vlmeta", "data"],
        default="tree",
        help="Panel to focus on startup",
    )
    parser.add_argument(
        "--mouse",
        action="store_true",
        help="Capture the mouse for clicking and scrolling (disables the terminal's native text selection)",
    )
    parser.add_argument(
        "--max",
        dest="maximized",
        action="store_true",
        help="Maximize the focused panel on startup (same as pressing 'm')",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        urlpath, download_url, info_url = resolve_source(args.urlpath, args.download)
    except ValueError as exc:
        parser.error(str(exc))

    import blosc2

    if blosc2.IS_WASM:
        print(
            "b2view is an interactive terminal UI and is not supported in the "
            "Pyodide/WebAssembly build of blosc2:\nthere is no terminal driver "
            "(termios) available in this environment.\n"
            "Run b2view from a native (CPython) install instead.",
            file=sys.stderr,
        )
        return 1

    try:
        from blosc2.b2view.app import B2ViewApp
    except ImportError as exc:
        print(
            'b2view could not import its TUI dependencies. Install them with:\n\n    pip install "blosc2[tui]"\n',
            file=sys.stderr,
        )
        print(f"Original import error: {exc}", file=sys.stderr)
        return 2

    app = B2ViewApp(
        urlpath,
        start_path=args.path,
        start_panel=args.panel,
        start_maximized=args.maximized,
        preview_rows=args.preview_rows,
        preview_cols=args.preview_cols,
        download_url=download_url,
        info_url=info_url,
    )
    app.run(mouse=args.mouse)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
