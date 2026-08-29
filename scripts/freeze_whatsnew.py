#!/usr/bin/env python3
"""
Freeze the executable examples in a whatsnew file.

Old whatsnew notes are a historical record, but their ``.. ipython::`` blocks are
re-executed on every doc build against whatever pandas is current.  That is how
they end up illustrating something other than the release they describe, and it
is why they periodically break the (warnings-as-errors) doc build when a
deprecation they demonstrate is enforced.

This script rewrites those blocks as static ``.. code-block:: pycon`` sessions.
By default the code and its output are taken from the release the file describes:
the source is read back out of git at that release's tag and executed against
that release's pandas, fetched on the fly with ``pixi exec``.  That is not quite
the byte-for-byte shipped page -- the directive shared one IPython shell across the
whole build, so display state leaked between pages -- but it is what the release
itself produced, and nothing about it can rot again.

    python scripts/freeze_whatsnew.py doc/source/whatsnew/v2.0.0.rst
    python scripts/freeze_whatsnew.py --diff doc/source/whatsnew/v1.1.0.rst

conda-forge has no pandas older than 0.17.1, so files below that can only be
rendered against the pandas in the current environment::

    python scripts/freeze_whatsnew.py --current doc/source/whatsnew/v0.13.0.rst

which stops the page from executing but does not make it historically accurate.

A block the sandbox cannot run -- one that needs the network or a package that will
not resolve alongside that old pandas -- is refused rather than frozen around a
fabricated traceback.  ``--skip-unusable`` leaves those blocks executable instead,
so the rest of the file can still be frozen.  One block needs both passes in turn:
v0.20.0's compressed-URL example points at ``branch='master'`` in the released
text, which 404s since the rename, so it is frozen from the current source with
``--current`` after the era pass has done the other 53::

    python scripts/freeze_whatsnew.py --skip-unusable <the file>
    python scripts/freeze_whatsnew.py --current --skip-unusable <the file>
"""

from __future__ import annotations

import argparse
import ast
import difflib
import json
import os
from platform import machine
import re
import subprocess
import sys
import tempfile
from typing import (
    TYPE_CHECKING,
    Any,
)
import warnings

if TYPE_CHECKING:
    from collections.abc import Sequence

# One statement, comment or magic inside a directive.
Item = dict[str, Any]
# What the era interpreter hands back: per-block, per-item output.
Results = dict[str, Any]

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOC = os.path.join(REPO, "doc")
WHATSNEW = os.path.join(DOC, "source", "whatsnew")

# Tolerant of a space before the colons.  docutils accepts that spelling and
# released files use it; the current tree does not, and a pre-commit hook now
# forbids it, so it only ever turns up in the revisions read back out of git.
DIRECTIVE = re.compile(r"^(\s*)\.\.\s+ipython\s*::\s*(.*?)\s*$")
OPTION = re.compile(r"^\s*:(\w+):\s*(.*)$")
DECORATOR = re.compile(r"^@(\w+)")
MAGIC = re.compile(r"^\s*[%!]")
PROMPT_IN = re.compile(r"^\s*In \[\d+\]:\s?(.*)$")
PROMPT_CONTINUE = re.compile(r"^\s*\.{3,}:\s?(.*)$")
VERSION = re.compile(r"^v(\d+)\.(\d+)\.(\d+)\.(?:rst|txt)$")
# Failures that say something about the sandbox rather than about the release,
# matched on the exception's type rather than on its formatted text -- the text
# puts SyntaxError behind a File line and wraps long messages, so any pattern
# over it goes stale the moment the formatting changes.
UNUSABLE_TYPES = frozenset(
    {
        "ConnectionError",
        "ConnectionResetError",
        "FileNotFoundError",
        "HTTPError",
        "ImportError",
        "IndentationError",
        "ModuleNotFoundError",
        "SyntaxError",
        "TimeoutError",
        "URLError",
        "gaierror",
        # Synthesised by freeze() rather than raised: see the raises loop.
        "Unparsed",
        "Volatile",
    }
)

# pandas 0.17-0.23 run on CPython 3.5/3.6, where a dict does not keep insertion
# order and its iteration order is salted per process.  An example that aggregates
# with a dict literal would otherwise freeze to whichever column order that run
# happened to produce, and the two-run check below would only catch it half the
# time.  Pinning the seed makes the output reproducible instead.
STABLE_ENV = {**os.environ, "PYTHONHASHSEED": "0"}

# conda-forge's oldest solvable pandas.
OLDEST_ERA = (0, 17)
# osx-arm64 packages only start at 1.1; older ones run under Rosetta.
OLDEST_ARM64 = (1, 1)

# doc/source/conf.py injects this into every page as {{ header }}.
HEADER = [
    "import numpy as np",
    "import pandas as pd",
    "np.random.seed(123456)",
    "np.set_printoptions(precision=4, suppress=True)",
    "pd.options.display.max_rows = 15",
]

# Older trees had no {{ header }}: whatsnew was one page that included every
# release's .txt, and its setup block lived at the top of that page.
LEGACY_HEADER_PAGE = "doc/source/whatsnew.rst"

# Optional dependencies an example may need, keyed by a marker in its source.
# Everything is fetched into the throwaway era environment, so this only has to
# be good enough to avoid an ImportError standing in for the real output.
OPTIONAL_DEPS = {
    "pyarrow": ["pyarrow"],
    "read_xml": ["lxml"],
    "to_xml": ["lxml"],
    "read_html": ["lxml", "html5lib", "beautifulsoup4"],
    ".style": ["jinja2"],
    "Styler": ["jinja2"],
    "matplotlib": ["matplotlib"],
    ".plot(": ["matplotlib"],
    "scipy": ["scipy"],
    "interpolate": ["scipy"],
    "read_excel": ["openpyxl", "xlrd"],
    "to_excel": ["openpyxl"],
    "HDFStore": ["pytables"],
    "_hdf": ["pytables"],
    "numexpr": ["numexpr"],
    "read_sql": ["sqlalchemy"],
    "to_sql": ["sqlalchemy"],
    "read_feather": ["pyarrow"],
    "read_parquet": ["pyarrow"],
    "to_parquet": ["pyarrow"],
    "to_markdown": ["tabulate"],
    "numba": ["numba"],
}

# Runs under the era interpreter, which may be as old as CPython 3.5.
# Keep it free of f-strings and any other modern syntax.
RUNNER = """
import io
import json
import sys
import traceback
import warnings


def main():
    payload = json.load(open(sys.argv[1]))
    namespace = {}
    for line in payload["header"]:
        exec(compile(line, "<header>", "exec"), namespace)

    results = []
    for block in payload["blocks"]:
        rendered = []
        for item in block:
            if item["kind"] != "code":
                rendered.append({"output": "", "error": None, "error_type": None})
                continue
            source = item["source"]
            buffer = io.StringIO()
            saved = sys.stdout
            sys.stdout = buffer
            error = None
            error_type = None
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                try:
                    try:
                        code = compile(source, "<doc>", "eval")
                    except SyntaxError:
                        exec(compile(source, "<doc>", "exec"), namespace)
                    else:
                        value = eval(code, namespace)
                        if value is not None:
                            namespace["_"] = value
                            buffer.write(repr(value) + "\\n")
                except BaseException:
                    kind, value = sys.exc_info()[:2]
                    error_type = kind.__name__
                    error = "".join(
                        traceback.format_exception_only(kind, value)
                    ).strip()
                    if error_type == "SyntaxError":
                        # format_exception_only prefixes the offending File line,
                        # which is a path inside the sandbox.
                        error = error.split("\\n")[-1]
                finally:
                    sys.stdout = saved
            rendered.append(
                {
                    "output": buffer.getvalue(),
                    "error": error,
                    "error_type": error_type,
                }
            )
        results.append(rendered)

    import pandas
    json.dump(
        {"results": results, "version": pandas.__version__},
        open(payload["outfile"], "w"),
    )


main()
"""


class Block:
    """One ``.. ipython::`` directive."""

    def __init__(
        self,
        start: int,
        stop: int,
        base: int,
        args: list[str],
        options: dict[str, str],
        content: list[str],
    ) -> None:
        self.start = start
        self.stop = stop
        self.base = base
        self.args = args
        self.options = options
        self.content = content

    @property
    def suppressed(self) -> bool:
        return "suppress" in self.options

    @property
    def verbatim(self) -> bool:
        return "verbatim" in self.options

    @property
    def key(self) -> str:
        """Identity for matching across revisions.

        Blackening the docs reflowed and requoted most of these examples, so
        neither layout nor quote style can be part of a block's identity.
        """
        return re.sub(r"\s+", "", " ".join(self.content).replace("'", '"'))


def parse_blocks(text: str) -> list[Block]:
    """Return every ``.. ipython::`` directive in ``text``, in source order."""
    lines = text.split("\n")
    blocks = []
    index = 0
    while index < len(lines):
        match = DIRECTIVE.match(lines[index])
        if match is None:
            index += 1
            continue
        base = len(match.group(1))
        args = match.group(2).split()
        stop = index + 1
        body = []
        while stop < len(lines):
            line = lines[stop]
            if not line.strip():
                body.append("")
                stop += 1
                continue
            if len(line) - len(line.lstrip()) <= base:
                break
            body.append(line)
            stop += 1
        while body and not body[-1].strip():
            body.pop()
            stop -= 1

        # Options come off first: they are indented relative to the directive,
        # not to the code, and several blocks indent them one column shallower.
        # Letting them into the dedent leaves every code line with a stray space,
        # which turns the whole body into unparsable text.
        options = {}
        while body:
            option = OPTION.match(body[0])
            if option is None:
                break
            options[option.group(1)] = option.group(2)
            body.pop(0)
        while body and not body[0].strip():
            body.pop(0)

        filled = [line for line in body if line.strip()]
        inner = min(len(ln) - len(ln.lstrip()) for ln in filled) if filled else base + 3
        body = [line[inner:] if line.strip() else "" for line in body]

        blocks.append(Block(index, stop, base, args, options, body))
        index = stop
    return blocks


def strip_prompts(content: list[str]) -> list[str]:
    """Undo an already-prompted block, the way the ipython directive does.

    A handful of blocks are written as ``In [5]: ...`` transcripts inside a live
    directive.  The directive re-executes the prompted lines and throws the
    pasted output away, so that is what has to be executed here too -- taking
    the text at face value would try to run ``In [5]: ...`` as Python.
    """
    if not any(PROMPT_IN.match(line) for line in content):
        return list(content)
    stripped = []
    for line in content:
        match = PROMPT_IN.match(line) or PROMPT_CONTINUE.match(line)
        if match is not None:
            stripped.append(match.group(1))
        elif not line.strip():
            stripped.append("")
    return stripped


def split_items(content: list[str]) -> list[Item]:
    """Split a directive body into statements, comments and magics.

    ``blank_before`` records whether the source separated this item from the
    previous one, so the frozen block can reproduce the author's spacing instead
    of the blank line the ipython directive inserts after every prompt.
    """
    lines = strip_prompts(content)
    decorators: dict[int, list[str]] = {}
    for number, line in enumerate(lines):
        match = DECORATOR.match(line.strip())
        if match is not None:
            decorators.setdefault(number, []).append(match.group(1))

    # ``@suppress`` and friends are not Python, and neither are ``%magics``.
    # Swap in placeholders so the block still parses, then map back by line.
    parseable = []
    magics = set()
    for number, line in enumerate(lines):
        if number in decorators:
            parseable.append("pass")
        elif MAGIC.match(line):
            magics.add(number)
            parseable.append("pass")
        else:
            parseable.append(line)

    try:
        with warnings.catch_warnings():
            # Released examples predate the escape-sequence warnings.
            warnings.simplefilter("ignore", SyntaxWarning)
            tree = ast.parse("\n".join(parseable))
    except SyntaxError:
        return [{"kind": "unparsed", "source": "\n".join(lines), "blank_before": False}]

    items: list[Item] = []
    pending: list[str] = []
    consumed = 0
    for node in tree.body:
        start = node.lineno - 1
        stop = node.end_lineno or start + 1
        if start in decorators:
            pending.extend(decorators[start])
            consumed = stop
            continue

        gap = lines[consumed:start]
        blank_before = any(not line.strip() for line in gap)
        for line in gap:
            if line.strip().startswith("#"):
                items.append(
                    {"kind": "comment", "source": line, "blank_before": blank_before}
                )
                blank_before = False

        kind = (
            "magic"
            if any(number in magics for number in range(start, stop))
            else "code"
        )
        statement = "\n".join(lines[start:stop])
        if kind == "code" and statement.rstrip().endswith(";"):
            # A trailing semicolon was IPython telling the directive not to show
            # the result -- usually because a screenshot followed.  It is not
            # Python, the screenshots are long gone, and leaving it in makes the
            # statement uncompilable as an expression, so the block would freeze
            # with no output at all.
            statement = statement.rstrip()[:-1]
        items.append(
            {
                "kind": kind,
                "source": statement,
                "blank_before": blank_before,
                "decorators": pending,
            }
        )
        pending = []
        consumed = stop

    items.extend(
        {"kind": "comment", "source": line, "blank_before": False}
        for line in lines[consumed:]
        if line.strip().startswith("#")
    )
    return items


def era_of(name: str) -> tuple[int, int]:
    """The ``(major, minor)`` release a whatsnew file describes."""
    match = VERSION.match(name)
    if match is None:
        raise ValueError(f"not a whatsnew file: {name}")
    return int(match.group(1)), int(match.group(2))


def era_tag(era: tuple[int, int]) -> str | None:
    """The last tag in ``era``'s series -- the file's final released state."""
    tags = subprocess.run(
        ["git", "-C", REPO, "tag"], capture_output=True, text=True, check=True
    ).stdout.split()
    pattern = re.compile(rf"^v{era[0]}\.{era[1]}\.(\d+)$")
    matched = [
        (int(match.group(1)), tag)
        for tag, match in ((tag, pattern.match(tag)) for tag in tags)
        if match is not None
    ]
    if not matched:
        return None
    return max(matched)[1]


def era_source(tag: str, name: str) -> str | None:
    """The whatsnew file as it stood at ``tag``.  It was ``.txt`` before 0.24."""
    stem = os.path.splitext(name)[0]
    for extension in (".rst", ".txt"):
        path = f"doc/source/whatsnew/{stem}{extension}"
        result = subprocess.run(
            ["git", "-C", REPO, "show", f"{tag}:{path}"],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode == 0:
            return result.stdout
    return None


SIMILAR_ENOUGH = 0.6


def era_header(tag: str, workdir: str) -> list[str]:
    """The setup the doc build ran before this release's whatsnew page.

    Getting this from the release matters: until 0.24 the whatsnew page did
    ``from pandas import *``, so its examples say ``Timestamp`` where today's
    say ``pd.Timestamp``.
    """
    lines: list[str] = []
    result = subprocess.run(
        ["git", "-C", REPO, "show", f"{tag}:{LEGACY_HEADER_PAGE}"],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode == 0:
        for block in parse_blocks(result.stdout):
            if block.suppressed:
                lines.extend(line for line in block.content if line.strip())
    if not lines:
        lines = list(HEADER)

    # Determinism and the data files the examples read, neither of which the
    # page's own setup provides.
    return [
        *lines,
        "import numpy as np",
        "np.random.seed(123456)",
        "np.set_printoptions(precision=4, suppress=True)",
        "import os",
        f"os.chdir({os.path.join(workdir, 'doc')!r})",
    ]


def era_doc_tree(tag: str) -> str:
    """Extract the release's ``doc/data`` so ``read_csv('data/...')`` still works."""
    workdir = tempfile.mkdtemp(prefix="freeze-whatsnew-doc-")
    archive = subprocess.run(
        ["git", "-C", REPO, "archive", tag, "doc/data"],
        capture_output=True,
        check=False,
    )
    if archive.returncode == 0:
        subprocess.run(["tar", "-x", "-C", workdir], input=archive.stdout, check=False)
    os.makedirs(os.path.join(workdir, "doc"), exist_ok=True)
    return workdir


def match_blocks(era: list[Block], current: list[Block]) -> dict[int, int]:
    """Map each current block onto its counterpart in the era source.

    The two rarely line up one-to-one: blocks frozen by earlier passes are still
    ``.. ipython::`` in the older revision, so the era file usually has more, and
    a decade of editing means the survivors are not identical either.  Align them
    the way a diff would -- a monotone pairing that maximises total similarity --
    rather than pairing by position.
    """
    era_keys = [block.key for block in era]
    current_keys = [block.key for block in current]

    def similarity(era_index: int, current_index: int) -> float:
        ratio = difflib.SequenceMatcher(
            None, era_keys[era_index], current_keys[current_index], autojunk=False
        ).ratio()
        return ratio if ratio >= SIMILAR_ENOUGH else 0.0

    rows, cols = len(era), len(current)
    best = [[0.0] * (cols + 1) for _ in range(rows + 1)]
    for row in range(1, rows + 1):
        for col in range(1, cols + 1):
            paired = similarity(row - 1, col - 1)
            best[row][col] = max(
                best[row - 1][col],
                best[row][col - 1],
                best[row - 1][col - 1] + paired if paired else 0.0,
            )

    mapping = {}
    row, col = rows, cols
    while row > 0 and col > 0:
        paired = similarity(row - 1, col - 1)
        if paired and best[row][col] == best[row - 1][col - 1] + paired:
            mapping[col - 1] = row - 1
            row, col = row - 1, col - 1
        elif best[row][col] == best[row - 1][col]:
            row -= 1
        else:
            col -= 1
    return mapping


def dependencies(blocks: list[Block]) -> list[str]:
    needed = set()
    for block in blocks:
        text = "\n".join(block.content)
        for marker, packages in OPTIONAL_DEPS.items():
            if marker in text:
                needed.update(packages)
    return sorted(needed)


def run_era(
    blocks: list[Block],
    prelude: list[Block],
    era: tuple[int, int],
    tag: str,
    extra: list[str],
    platform: str | None,
    verbose: bool,
) -> Results:
    """Execute ``blocks`` against ``era``'s pandas, fetched with ``pixi exec``."""
    payload = {
        "header": era_header(tag, era_doc_tree(tag)),
        "blocks": [split_items(block.content) for block in prelude + blocks],
    }
    results = _execute(payload, era, extra, platform, verbose)
    results["results"] = results["results"][len(prelude) :]
    return results


def mark_volatile(first: Results, second: Results) -> Results:
    """Flag every item whose output is not the same twice running.

    Enumerating the ways output can be volatile does not work -- an object's
    address looks nothing like a wall clock, which looks nothing like an ``id()``
    used as a label -- so the file is simply run twice, from scratch each time,
    and anything that disagrees with itself is refused rather than frozen.
    """
    for left, right in zip(first["results"], second["results"], strict=True):
        for item, other in zip(left, right, strict=True):
            item["volatile"] = item["output"] != other["output"]
    return first


def _execute(
    payload: dict[str, Any],
    era: tuple[int, int],
    extra: list[str],
    platform: str | None,
    verbose: bool,
) -> Results:
    directory = tempfile.mkdtemp(prefix="freeze-whatsnew-")
    runner = os.path.join(directory, "runner.py")
    inputs = os.path.join(directory, "payload.json")
    outputs = os.path.join(directory, "results.json")
    payload["outfile"] = outputs
    with open(runner, "w") as handle:
        handle.write(RUNNER)
    with open(inputs, "w") as handle:
        json.dump(payload, handle)

    spec = ["pandas=={}.{}.*".format(*era)]
    # conda-forge's pandas metadata does not cap numpy, so an unpinned solve
    # pairs old pandas with a numpy that breaks it: >=1.24 removed np.bool,
    # which these versions use at import, and >=1.20 breaks their compiled
    # dtype handling ("Cannot interpret ... as a data type").
    if era < (1, 0):
        spec.append("numpy<1.20")
    elif era < (1, 5):
        spec.append("numpy<1.24")
    spec.extend(extra)

    command = ["pixi", "exec"]
    if platform:
        command += ["--platform", platform]
    for item in spec:
        command += ["-s", item]
    command += ["python", runner, inputs]

    if verbose:
        print("  $ {}".format(" ".join(command)), file=sys.stderr)

    def once() -> Results:
        result = subprocess.run(
            command,
            capture_output=not verbose,
            text=True,
            check=False,
            env=STABLE_ENV,
        )
        if result.returncode != 0:
            detail = "" if verbose else f"\n{result.stderr}"
            raise RuntimeError(
                "pixi exec failed for pandas {}.{}{}".format(*era, detail)
            )
        with open(outputs) as handle:
            return json.load(handle)

    return mark_volatile(once(), once())


def run_current(blocks: list[Block]) -> Results:
    """Execute ``blocks`` in this interpreter, for files with no era pandas."""
    payload = {
        "header": [*HEADER, "import os", f"os.chdir({DOC!r})"],
        "blocks": [split_items(block.content) for block in blocks],
    }
    directory = tempfile.mkdtemp(prefix="freeze-whatsnew-")
    runner = os.path.join(directory, "runner.py")
    inputs = os.path.join(directory, "payload.json")
    outputs = os.path.join(directory, "results.json")
    payload["outfile"] = outputs
    with open(runner, "w") as handle:
        handle.write(RUNNER)
    with open(inputs, "w") as handle:
        json.dump(payload, handle)

    def once() -> Results:
        subprocess.run([sys.executable, runner, inputs], check=True, env=STABLE_ENV)
        with open(outputs) as handle:
            return json.load(handle)

    return mark_volatile(once(), once())


def normalize_backticks(comment: str) -> str:
    """Double up single-backtick spans in a restored comment.

    The rst-backticks hook cannot see that a line is inside a literal block, and
    several of these comments were double-backticked years after release to keep
    it quiet.  Re-applying that here keeps the release text honest -- it is a
    comment, not something pandas ever printed -- without reintroducing a lint
    failure the tree has already fixed once.
    """
    return re.sub(r"(?<!`)`([^`\n]+)`(?!`)", r"``\1``", comment)


def render(items: list[Item], results: list[dict[str, Any]]) -> list[str]:
    """Turn one executed block into doctest-style lines."""
    out: list[str] = []
    for item, result in zip(items, results, strict=True):
        if item["blank_before"] and out:
            out.append("")

        if item["kind"] == "comment":
            out.extend(normalize_backticks(item["source"]).split("\n"))
            continue

        if "suppress" in item.get("decorators", []):
            # The directive hid the input too; it ran only to set something up.
            continue

        source = item["source"].split("\n")
        out.append(">>> " + source[0])
        out.extend(("... " + line) if line else "..." for line in source[1:])

        if item["kind"] == "magic":
            continue
        if result["error"]:
            out.append("Traceback (most recent call last):")
            out.append("    ...")
            out.extend(result["error"].split("\n"))
        elif result["output"]:
            out.extend(result["output"].rstrip("\n").split("\n"))
    return out


def freeze(
    path: str,
    current_pandas: bool,
    platform: str | None,
    skip_unusable: bool,
    verbose: bool,
) -> tuple[str | None, dict[str, Any]]:
    """Return the rewritten text for ``path``, plus a report on what happened."""
    name = os.path.basename(path)
    with open(path) as handle:
        text = handle.read()
    blocks = parse_blocks(text)
    if not blocks:
        return None, {"blocks": 0, "reason": "nothing to freeze"}

    if all(not block.content for block in blocks):
        return None, {"blocks": len(blocks), "reason": "no directive has a body"}

    era = era_of(name)
    report: dict[str, Any] = {
        "blocks": len(blocks),
        "era": era,
        "restored": 0,
        "unmatched": 0,
    }
    restored = 0

    if current_pandas or era < OLDEST_ERA:
        if not current_pandas:
            report["reason"] = "no pandas {}.{} on conda-forge".format(*era)
            return None, report
        sources = blocks
        mapping = {index: index for index in range(len(blocks))}
        results = run_current(blocks)
        report["rendered_with"] = results["version"] + " (current)"
    else:
        tag = era_tag(era)
        original = era_source(tag, name) if tag is not None else None
        if tag is None or original is None:
            report["reason"] = f"no released copy of {name} to read back"
            return None, report
        sources = parse_blocks(original)
        mapping = match_blocks(sources, blocks)
        # Suppressed blocks are dropped rather than rendered, so a setup block
        # added to the page after the release need not exist in the era source.
        orphans = [
            index
            for index, block in enumerate(blocks)
            if index not in mapping and not block.suppressed
        ]
        # A suppressed block with no era counterpart is a shim added later to
        # keep the page runnable -- most often ``from pandas import *``, which
        # the doc header used to do for everyone.  Run it first to put the era's
        # ambient namespace back; it renders nothing either way.
        prelude = [
            block
            for index, block in enumerate(blocks)
            if index not in mapping and block.suppressed
        ]
        report["unmatched"] = len(orphans)
        if orphans:
            report["reason"] = "{} block(s) have no counterpart at {} ({})".format(
                len(orphans),
                tag,
                ", ".join(f"line {blocks[index].start + 1}" for index in orphans),
            )
            return None, report
        if platform is None and era < OLDEST_ARM64 and _is_arm64():
            platform = "osx-64"
        results = run_era(
            sources,
            prelude,
            era,
            tag,
            dependencies(sources + prelude),
            platform,
            verbose,
        )
        report["tag"] = tag
        report["rendered_with"] = results["version"]

    # An example that raises is usually the point of the block, but an example
    # that raises because the sandbox lacks a package or a network is not
    # something to bake into a release note.
    raises = []
    volatile: list[tuple[int, str]] = []
    # Blocks rendered as input with no output: an IPython magic, or a statement
    # whose output is not the same twice running.  Leaving such a block
    # executable is not an option -- once its setup blocks are frozen it no
    # longer has the state it needs, and the doc build fails on it.
    input_only = set()
    for index, block in enumerate(blocks):
        if block.suppressed or block.verbatim:
            continue
        items = split_items(sources[mapping[index]].content)
        for item, result in zip(items, results["results"][mapping[index]], strict=True):
            if item["kind"] == "magic":
                input_only.add(index)
            elif item["kind"] == "unparsed":
                # The body did not parse as Python, so nothing was executed.
                # Rendering it anyway would emit a promptless, output-less
                # dump that still looks like a frozen session.
                raises.append(
                    (
                        block.start + 1,
                        item["source"].split("\n")[0].strip(),
                        "block did not parse as Python",
                        "Unparsed",
                    )
                )
            elif result["error"]:
                raises.append(
                    (
                        block.start + 1,
                        item["source"].split("\n")[0],
                        result["error"],
                        result["error_type"],
                    )
                )
            elif result.get("volatile"):
                input_only.add(index)
                volatile.append((block.start + 1, item["source"].split("\n")[0]))
    report["volatile"] = volatile
    unusable = {entry[0] for entry in raises if entry[3] in UNUSABLE_TYPES}
    if unusable and not skip_unusable:
        first: dict[int, str] = {}
        for line, _, error, _kind in raises:
            first.setdefault(line, error)
        report["reason"] = (
            "the sandbox could not run {} block(s): {}; pass --skip-unusable to "
            "leave them executable".format(
                len(unusable),
                "; ".join(
                    f"line {line} -- {error}"
                    for line, error in sorted(first.items())
                    if line in unusable
                ),
            )
        )
        return None, report
    raises = [entry for entry in raises if entry[0] not in unusable]
    report["skipped"] = sorted(unusable)

    lines = text.split("\n")
    for index in reversed(range(len(blocks))):
        block = blocks[index]
        pad = " " * block.base

        if block.start + 1 in unusable:
            # Left executable on purpose; see --skip-unusable.
            continue

        if block.suppressed:
            # Setup and teardown for code that no longer runs.  Take one of the
            # blank lines that surrounded it too, so the paragraph either side
            # keeps exactly the spacing it had.
            stop = block.stop
            before_blank = block.start > 0 and not lines[block.start - 1].strip()
            after_blank = stop < len(lines) and not lines[stop].strip()
            if before_blank and after_blank:
                stop += 1
            filler = [] if before_blank or after_blank else [""]
            lines[block.start : stop] = filler
            continue

        source = sources[mapping[index]]
        if block.verbatim or index in input_only:
            # Hand-written transcripts; the ``%timeit`` pair in v0.13.0, which is
            # not Python and only ever printed a timing of whichever machine
            # built the docs; and blocks whose output is not the same twice
            # running (a wall clock, an ``id()``, a dict order under CPython
            # 3.5).  Show what to run rather than pin a value that was never
            # reproducible.
            body = [".. code-block:: ipython", ""]
            body += ["   " + line if line else "" for line in block.content]
        else:
            items = split_items(source.content)
            rendered = render(items, results["results"][mapping[index]])
            if not rendered:
                # Every statement was suppressed, or the era block was empty.
                # An empty literal block is a docutils error, so drop it the way
                # a wholly-suppressed directive is dropped.
                report.setdefault("emptied", []).append(block.start + 1)
                lines[block.start : block.stop] = []
                continue
            body = [".. code-block:: pycon", ""]
            body += ["   " + line if line else "" for line in rendered]
            if source.key != block.key:
                restored += 1

        replacement = [(pad + line).rstrip() if line else "" for line in body]
        if block.stop < len(lines) and lines[block.stop].strip():
            # The directive swallowed its own content, so a few blocks run
            # straight into the next paragraph.  A literal block may not.
            replacement.append("")
        lines[block.start : block.stop] = replacement

    report["restored"] = restored
    report["raises"] = raises
    return "\n".join(lines), report


def _is_arm64() -> bool:
    return machine() in ("arm64", "aarch64")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", help="whatsnew files to freeze")
    parser.add_argument(
        "--current",
        action="store_true",
        help="render against the installed pandas instead of the release's",
    )
    parser.add_argument(
        "--diff", action="store_true", help="print a diff instead of writing"
    )
    parser.add_argument("--platform", help="conda platform to fetch the era pandas for")
    parser.add_argument(
        "--skip-unusable",
        action="store_true",
        help="leave blocks the sandbox cannot run (network, missing package) "
        "executable instead of refusing the file",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    failed = False
    for path in args.paths:
        name = os.path.basename(path)
        try:
            frozen, report = freeze(
                path, args.current, args.platform, args.skip_unusable, args.verbose
            )
        except Exception as err:
            print(f"{name}: {err}", file=sys.stderr)
            failed = True
            continue

        if frozen is None:
            print("{}: skipped -- {}".format(name, report.get("reason", "?")))
            failed = True
            continue

        summary = "{}: {} block(s) rendered with pandas {}".format(
            name, report["blocks"], report["rendered_with"]
        )
        if report.get("tag"):
            summary += " from {}".format(report["tag"])
        if report["restored"]:
            summary += ", {} restored to their released text".format(report["restored"])
        if report["raises"]:
            summary += ", {} example(s) raise".format(len(report["raises"]))
        if report.get("skipped"):
            summary += ", {} block(s) left executable".format(len(report["skipped"]))
        if report.get("volatile"):
            summary += ", {} shown without output".format(len(report["volatile"]))
        print(summary)
        for line, statement in report.get("volatile", []):
            print(f"    line {line}: {statement} -> output differs between runs")
        for line in report.get("skipped", []):
            print(f"    line {line}: left as a live directive (--skip-unusable)")
        for line, statement, error, _kind in report["raises"]:
            print(f"    line {line}: {statement} -> {error}")

        if args.diff:
            with open(path) as handle:
                before = handle.read()
            sys.stdout.writelines(
                difflib.unified_diff(
                    before.splitlines(True),
                    frozen.splitlines(True),
                    fromfile=path,
                    tofile=path + " (frozen)",
                )
            )
        else:
            with open(path, "w") as handle:
                handle.write(frozen)

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
