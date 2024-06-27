# Licensed under the Apache License: http://www.apache.org/licenses/LICENSE-2.0
# For details: https://github.com/nedbat/coveragepy/blob/master/NOTICE.txt

"""Results of coverage measurement."""

from __future__ import annotations

import collections
import dataclasses

from collections.abc import Container
from typing import Iterable, TYPE_CHECKING

from coverage.exceptions import ConfigError
from coverage.misc import nice_pair
from coverage.types import TArc, TLineNo

if TYPE_CHECKING:
    from coverage.data import CoverageData
    from coverage.plugin import FileReporter


def analysis_from_file_reporter(
    data: CoverageData,
    precision: int,
    file_reporter: FileReporter,
    filename: str,
) -> Analysis:
    """Create an Analysis from a FileReporter."""
    has_arcs = data.has_arcs()
    statements = file_reporter.lines()
    excluded = file_reporter.excluded_lines()
    executed = file_reporter.translate_lines(data.lines(filename) or [])

    if has_arcs:
        _arc_possibilities_set = file_reporter.arcs()
        _arcs_executed_set = file_reporter.translate_arcs(data.arcs(filename) or [])
        exit_counts = file_reporter.exit_counts()
        no_branch = file_reporter.no_branch_lines()
    else:
        _arc_possibilities_set = set()
        _arcs_executed_set = set()
        exit_counts = {}
        no_branch = set()

    return Analysis(
        precision=precision,
        filename=filename,
        has_arcs=has_arcs,
        statements=statements,
        excluded=excluded,
        executed=executed,
        _arc_possibilities_set=_arc_possibilities_set,
        _arcs_executed_set=_arcs_executed_set,
        exit_counts=exit_counts,
        no_branch=no_branch,
    )


@dataclasses.dataclass
class Analysis:
    """The results of analyzing a FileReporter."""

    precision: int
    filename: str
    has_arcs: bool
    statements: set[TLineNo]
    excluded: set[TLineNo]
    executed: set[TLineNo]
    _arc_possibilities_set: set[TArc]
    _arcs_executed_set: set[TArc]
    exit_counts: dict[TLineNo, int]
    no_branch: set[TLineNo]

    def __post_init__(self) -> None:
        self.arc_possibilities = sorted(self._arc_possibilities_set)
        self.arcs_executed = sorted(self._arcs_executed_set)
        self.missing = self.statements - self.executed

        if self.has_arcs:
            n_branches = self._total_branches()
            mba = self.missing_branch_arcs()
            n_partial_branches = sum(len(v) for k,v in mba.items() if k not in self.missing)
            n_missing_branches = sum(len(v) for k,v in mba.items())
        else:
            n_branches = n_partial_branches = n_missing_branches = 0

        self.numbers = Numbers(
            precision=self.precision,
            n_files=1,
            n_statements=len(self.statements),
            n_excluded=len(self.excluded),
            n_missing=len(self.missing),
            n_branches=n_branches,
            n_partial_branches=n_partial_branches,
            n_missing_branches=n_missing_branches,
        )

    def narrow(self, lines: Container[TLineNo]) -> Analysis:
        """Create a narrowed Analysis.

        The current analysis is copied to make a new one that only considers
        the lines in `lines`.
        """

        statements = {lno for lno in self.statements if lno in lines}
        excluded = {lno for lno in self.excluded if lno in lines}
        executed = {lno for lno in self.executed if lno in lines}

        if self.has_arcs:
            _arc_possibilities_set = {
                (a, b) for a, b in self._arc_possibilities_set
                if a in lines or b in lines
            }
            _arcs_executed_set = {
                (a, b) for a, b in self._arcs_executed_set
                if a in lines or b in lines
            }
            exit_counts = {
                lno: num for lno, num in self.exit_counts.items()
                if lno in lines
            }
            no_branch = {lno for lno in self.no_branch if lno in lines}
        else:
            _arc_possibilities_set = set()
            _arcs_executed_set = set()
            exit_counts = {}
            no_branch = set()

        return Analysis(
            precision=self.precision,
            filename=self.filename,
            has_arcs=self.has_arcs,
            statements=statements,
            excluded=excluded,
            executed=executed,
            _arc_possibilities_set=_arc_possibilities_set,
            _arcs_executed_set=_arcs_executed_set,
            exit_counts=exit_counts,
            no_branch=no_branch,
        )

    def missing_formatted(self, branches: bool = False) -> str:
        """The missing line numbers, formatted nicely.

        Returns a string like "1-2, 5-11, 13-14".

        If `branches` is true, includes the missing branch arcs also.

        """
        if branches and self.has_arcs:
            arcs = self.missing_branch_arcs().items()
        else:
            arcs = None

        return format_lines(self.statements, self.missing, arcs=arcs)

    def arcs_missing(self) -> list[TArc]:
        """Returns a sorted list of the un-executed arcs in the code."""
        missing = (
            p for p in self.arc_possibilities
                if p not in self.arcs_executed
                    and p[0] not in self.no_branch
                    and p[1] not in self.excluded
        )
        return sorted(missing)

    def arcs_unpredicted(self) -> list[TArc]:
        """Returns a sorted list of the executed arcs missing from the code."""
        # Exclude arcs here which connect a line to itself.  They can occur
        # in executed data in some cases.  This is where they can cause
        # trouble, and here is where it's the least burden to remove them.
        # Also, generators can somehow cause arcs from "enter" to "exit", so
        # make sure we have at least one positive value.
        unpredicted = (
            e for e in self.arcs_executed
                if e not in self.arc_possibilities
                    and e[0] != e[1]
                    and (e[0] > 0 or e[1] > 0)
        )
        return sorted(unpredicted)

    def _branch_lines(self) -> list[TLineNo]:
        """Returns a list of line numbers that have more than one exit."""
        return [l1 for l1,count in self.exit_counts.items() if count > 1]

    def _total_branches(self) -> int:
        """How many total branches are there?"""
        return sum(count for count in self.exit_counts.values() if count > 1)

    def missing_branch_arcs(self) -> dict[TLineNo, list[TLineNo]]:
        """Return arcs that weren't executed from branch lines.

        Returns {l1:[l2a,l2b,...], ...}

        """
        missing = self.arcs_missing()
        branch_lines = set(self._branch_lines())
        mba = collections.defaultdict(list)
        for l1, l2 in missing:
            if l1 in branch_lines:
                mba[l1].append(l2)
        return mba

    def executed_branch_arcs(self) -> dict[TLineNo, list[TLineNo]]:
        """Return arcs that were executed from branch lines.

        Returns {l1:[l2a,l2b,...], ...}

        """
        branch_lines = set(self._branch_lines())
        eba = collections.defaultdict(list)
        for l1, l2 in self.arcs_executed:
            if l1 in branch_lines:
                eba[l1].append(l2)
        return eba

    def branch_stats(self) -> dict[TLineNo, tuple[int, int]]:
        """Get stats about branches.

        Returns a dict mapping line numbers to a tuple:
        (total_exits, taken_exits).
        """

        missing_arcs = self.missing_branch_arcs()
        stats = {}
        for lnum in self._branch_lines():
            exits = self.exit_counts[lnum]
            missing = len(missing_arcs[lnum])
            stats[lnum] = (exits, exits - missing)
        return stats


@dataclasses.dataclass
class Numbers:
    """The numerical results of measuring coverage.

    This holds the basic statistics from `Analysis`, and is used to roll
    up statistics across files.

    """

    precision: int = 0
    n_files: int = 0
    n_statements: int = 0
    n_excluded: int = 0
    n_missing: int = 0
    n_branches: int = 0
    n_partial_branches: int = 0
    n_missing_branches: int = 0

    @property
    def n_executed(self) -> int:
        """Returns the number of executed statements."""
        return self.n_statements - self.n_missing

    @property
    def n_executed_branches(self) -> int:
        """Returns the number of executed branches."""
        return self.n_branches - self.n_missing_branches

    @property
    def pc_covered(self) -> float:
        """Returns a single percentage value for coverage."""
        if self.n_statements > 0:
            numerator, denominator = self.ratio_covered
            pc_cov = (100.0 * numerator) / denominator
        else:
            pc_cov = 100.0
        return pc_cov

    @property
    def pc_covered_str(self) -> str:
        """Returns the percent covered, as a string, without a percent sign.

        Note that "0" is only returned when the value is truly zero, and "100"
        is only returned when the value is truly 100.  Rounding can never
        result in either "0" or "100".

        """
        return display_covered(self.pc_covered, self.precision)

    @property
    def ratio_covered(self) -> tuple[int, int]:
        """Return a numerator and denominator for the coverage ratio."""
        numerator = self.n_executed + self.n_executed_branches
        denominator = self.n_statements + self.n_branches
        return numerator, denominator

    def __add__(self, other: Numbers) -> Numbers:
        return Numbers(
            self.precision,
            self.n_files + other.n_files,
            self.n_statements + other.n_statements,
            self.n_excluded + other.n_excluded,
            self.n_missing + other.n_missing,
            self.n_branches + other.n_branches,
            self.n_partial_branches + other.n_partial_branches,
            self.n_missing_branches + other.n_missing_branches,
        )

    def __radd__(self, other: int) -> Numbers:
        # Implementing 0+Numbers allows us to sum() a list of Numbers.
        assert other == 0   # we only ever call it this way.
        return self


def display_covered(pc: float, precision: int) -> str:
    """Return a displayable total percentage, as a string.

    Note that "0" is only returned when the value is truly zero, and "100"
    is only returned when the value is truly 100.  Rounding can never
    result in either "0" or "100".

    """
    near0 = 1.0 / 10 ** precision
    if 0 < pc < near0:
        pc = near0
    elif (100.0 - near0) < pc < 100:
        pc = 100.0 - near0
    else:
        pc = round(pc, precision)
    return "%.*f" % (precision, pc)


def _line_ranges(
    statements: Iterable[TLineNo],
    lines: Iterable[TLineNo],
) -> list[tuple[TLineNo, TLineNo]]:
    """Produce a list of ranges for `format_lines`."""
    statements = sorted(statements)
    lines = sorted(lines)

    pairs = []
    start = None
    lidx = 0
    for stmt in statements:
        if lidx >= len(lines):
            break
        if stmt == lines[lidx]:
            lidx += 1
            if not start:
                start = stmt
            end = stmt
        elif start:
            pairs.append((start, end))
            start = None
    if start:
        pairs.append((start, end))
    return pairs


def format_lines(
    statements: Iterable[TLineNo],
    lines: Iterable[TLineNo],
    arcs: Iterable[tuple[TLineNo, list[TLineNo]]] | None = None,
) -> str:
    """Nicely format a list of line numbers.

    Format a list of line numbers for printing by coalescing groups of lines as
    long as the lines represent consecutive statements.  This will coalesce
    even if there are gaps between statements.

    For example, if `statements` is [1,2,3,4,5,10,11,12,13,14] and
    `lines` is [1,2,5,10,11,13,14] then the result will be "1-2, 5-11, 13-14".

    Both `lines` and `statements` can be any iterable. All of the elements of
    `lines` must be in `statements`, and all of the values must be positive
    integers.

    If `arcs` is provided, they are (start,[end,end,end]) pairs that will be
    included in the output as long as start isn't in `lines`.

    """
    line_items = [(pair[0], nice_pair(pair)) for pair in _line_ranges(statements, lines)]
    if arcs is not None:
        line_exits = sorted(arcs)
        for line, exits in line_exits:
            for ex in sorted(exits):
                if line not in lines and ex not in lines:
                    dest = (ex if ex > 0 else "exit")
                    line_items.append((line, f"{line}->{dest}"))

    ret = ", ".join(t[-1] for t in sorted(line_items))
    return ret


def should_fail_under(total: float, fail_under: float, precision: int) -> bool:
    """Determine if a total should fail due to fail-under.

    `total` is a float, the coverage measurement total. `fail_under` is the
    fail_under setting to compare with. `precision` is the number of digits
    to consider after the decimal point.

    Returns True if the total should fail.

    """
    # We can never achieve higher than 100% coverage, or less than zero.
    if not (0 <= fail_under <= 100.0):
        msg = f"fail_under={fail_under} is invalid. Must be between 0 and 100."
        raise ConfigError(msg)

    # Special case for fail_under=100, it must really be 100.
    if fail_under == 100.0 and total != 100.0:
        return True

    return round(total, precision) < fail_under
