"""Run unit tests in parallel threads."""

# Based on unittest-ft, which is copyright Amethyst Reese, MIT license.

from __future__ import annotations

import os
import sys
import time
import random
import logging
import argparse
import unittest
import threading
import collections
from typing import Any, Self, TextIO
from dataclasses import dataclass, field
from collections.abc import Generator
from concurrent.futures import (
    Future,
    ThreadPoolExecutor,
    as_completed,
)

LOG = logging.getLogger(__name__)

DEFAULT_THREADS = os.cpu_count() or 4


class FTTestResult(unittest.TestResult):
    def __init__(
        self,
        stream: TextIO | None = None,
        descriptions: bool | None = None,
        verbosity: int | None = None,
    ) -> None:
        super().__init__(
            stream=stream, descriptions=descriptions, verbosity=verbosity
        )
        self.verbosity = verbosity or 1
        self.before = time.monotonic_ns()
        self.duration = 0
        self.collected_duration = 0

    def stopTest(self, test: Any) -> None:
        super().stopTest(test)
        self.duration = time.monotonic_ns() - self.before

    def stopTestRun(self) -> None:
        super().stopTestRun()
        self.duration = time.monotonic_ns() - self.before

    def __str__(self) -> str:
        items: dict[tuple[str, str], int] = collections.defaultdict(int)
        for test_case, trace in self.errors:
            items[f"ERROR: {test_case}", trace] += 1
        for test_case, trace in self.failures:
            items[f"FAIL: {test_case}", trace] += 1

        results = {
            f"{label}": trace for (label, trace), count in items.items()
        }
        longest = max(len(label) for label in results) if results else 70

        msg = "\n"
        msg += "\n".join(
            f"{'=' * longest}\n{label}\n{'-' * longest}\n{trace}"
            for label, trace in results.items()
        )
        msg += "-" * longest
        msg += f"\nRan {self.testsRun} tests in {format_ns(self.duration)}"

        saved = self.collected_duration - self.duration
        if saved > 0 and (saved / self.duration) > 0.10:
            msg += f" (saved {format_ns(self.collected_duration - self.duration)})"
        msg += "\n\n"

        msg += "OK" if self.wasSuccessful() else "FAILED"

        parts = []
        if self.errors:
            parts += [f"errors={len(self.errors)}"]
        if self.failures:
            parts += [f"failures={len(self.failures)}"]
        if self.skipped:
            parts += [f"skipped={len(self.skipped)}"]
        if self.expectedFailures:
            parts += [f"expected failures={len(self.expectedFailures)}"]

        if parts:
            msg += f" ({', '.join(parts)})"

        return msg

    def __add__(self, other: object) -> FTTestResult:
        if not isinstance(other, unittest.TestResult):
            return NotImplemented
        result = FTTestResult()
        result.errors = self.errors + other.errors
        result.expectedFailures = (
            self.expectedFailures + other.expectedFailures
        )
        result.failures = self.failures + other.failures
        result.skipped = self.skipped + other.skipped
        result.testsRun = self.testsRun + other.testsRun
        result.unexpectedSuccesses = (
            self.unexpectedSuccesses + other.unexpectedSuccesses
        )
        if isinstance(other, FTTestResult):
            result.collected_duration = self.duration + other.duration
        return result

    def __iadd__(self, other: object) -> Self:
        if not isinstance(other, unittest.TestResult):
            return NotImplemented
        self.errors += other.errors
        self.expectedFailures += other.expectedFailures
        self.failures += other.failures
        self.skipped += other.skipped
        self.testsRun += other.testsRun
        self.unexpectedSuccesses += other.unexpectedSuccesses
        if isinstance(other, FTTestResult):
            self.collected_duration += other.duration
        return self


def get_individual_tests(
    suite: unittest.TestSuite,
) -> Generator[unittest.TestCase]:
    for test in suite:
        if isinstance(test, unittest.TestSuite):
            yield from get_individual_tests(test)
        else:
            yield test


def run_single_test(suite: unittest.TestSuite) -> tuple[str, FTTestResult]:
    test_id = suite.id()
    LOG.debug(f"Running test {threading.get_ident()} {test_id}")
    result = FTTestResult(descriptions=True, verbosity=2)
    suite.run(result)
    LOG.debug("Finished test %s", test_id)
    return (test_id, result)


def format_ns(duration: int) -> str:
    if duration < 1_000_000_000:
        return f"{duration / 1_000_000:.2f}ms"
    else:
        return f"{duration / 1_000_000_000:.3f}s"


@dataclass
class Output:
    total: int
    futures: dict[Future[tuple[str, FTTestResult]], str] = field(
        default_factory=dict
    )
    stream: TextIO = sys.stdout
    verbosity: int = 1

    def __post_init__(self) -> None:
        self.count = 0

    def render(self, test_id: str, test_result: FTTestResult) -> None:
        stream = self.stream
        verbosity = self.verbosity

        self.count += 1
        if verbosity == 2:
            stream.write(
                f"[{self.count}/{self.total}] {test_id}"
                f" ... {'OK' if test_result.wasSuccessful() else 'FAIL'} "
                f" {format_ns(test_result.duration)}\n"
            )
        elif verbosity == 1:
            if test_result.errors:
                stream.write("E")
            elif test_result.failures:
                stream.write("F")
            elif test_result.expectedFailures:
                stream.write("x")
            elif test_result.skipped:
                stream.write("s")
            else:
                stream.write(".")
        stream.flush()


_EXCLUDE_CASES = set("""
    tables.filters.Filters
    tables.misc.enum
    tables.tests.test_array.SI1NACloseTestCase
    tables.tests.test_array.SI1NAOpenTestCase
    tables.tests.test_array.SI2NACloseTestCase
    tables.tests.test_array.SI2NAOpenTestCase
    tables.tests.test_basics.HDF5ErrorHandling
    tables.tests.test_basics.OpenFileFailureTestCase
    tables.tests.test_create.SetBloscMaxThreadsTestCase
    tables.tests.test_tablesMD.CompressTwoTablesTestCase
    tables.tests.test_utils.ptdumpTestCase.test_paths_windows
    tables.tests.test_utils.ptrepackTestCase.test_paths_windows
    tables.tests.test_utils.pttreeTestCase.test_paths_windows
    """.strip().split())


def _match_case(test_id: str) -> bool:
    for pat in _EXCLUDE_CASES:
        if test_id.startswith(pat):
            return True
    return False


def run_suite(
    suite: unittest.TestSuite,
    *,
    randomize: bool = False,
    threads: int = DEFAULT_THREADS,
    verbosity: int = 1,
    max_tests: int = 0,
) -> unittest.TestResult:
    test_cases = {}
    for test_case in get_individual_tests(suite):
        test_id = test_case.id()
        if _match_case(test_id):
            continue
        test_cases[test_id] = test_case
    test_ids = list(test_cases)
    if randomize:
        rnd = random.SystemRandom()
        rnd.shuffle(test_ids)
    else:
        test_ids.sort()
    if max_tests:
        test_ids = test_ids[:max_tests]

    LOG.info(
        "Ready to run %d tests:\n  %s", len(test_ids), "\n  ".join(test_ids)
    )

    output = Output(total=len(test_ids), verbosity=verbosity)
    result = FTTestResult()

    with ThreadPoolExecutor(max_workers=threads) as pool:
        futures = [
            pool.submit(run_single_test, test_cases[test_id])
            for test_id in test_ids
        ]
        for fut in as_completed(futures):
            test_id, test_result = fut.result()
            result += test_result
            output.render(test_id, test_result)
        result.stopTestRun()

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Run PyTables tests with free-threading support"
    )
    parser.add_argument(
        "--randomize",
        action="store_true",
        help="Run tests in random order",
    )
    parser.add_argument(
        "--verbosity",
        type=int,
        default=1,
        help="Logger verbosity level (default: 1)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=DEFAULT_THREADS,
        help=f"Number of worker threads to use (default: {DEFAULT_THREADS})",
    )
    parser.add_argument(
        "--max-tests",
        type=int,
        default=0,
        help="Maximum number of tests to run, 0 = unlimited (default: 0)",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbosity > 1 else logging.INFO
    )

    from tables.tests import test_suite

    suite = test_suite.suite()
    result = run_suite(
        suite,
        verbosity=args.verbosity,
        threads=args.threads,
        randomize=args.randomize,
        max_tests=args.max_tests,
    )
    print(result)


if __name__ == "__main__":
    main()
