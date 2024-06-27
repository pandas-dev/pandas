"""Plugin built-in to Flake8 to treat pyflakes as a plugin."""
from __future__ import annotations

import argparse
import ast
import logging
import os
from typing import Any
from typing import Generator

import pyflakes.checker

from flake8 import utils
from flake8.options.manager import OptionManager

LOG = logging.getLogger(__name__)

FLAKE8_PYFLAKES_CODES = {
    "UnusedImport": "F401",
    "ImportShadowedByLoopVar": "F402",
    "ImportStarUsed": "F403",
    "LateFutureImport": "F404",
    "ImportStarUsage": "F405",
    "ImportStarNotPermitted": "F406",
    "FutureFeatureNotDefined": "F407",
    "PercentFormatInvalidFormat": "F501",
    "PercentFormatExpectedMapping": "F502",
    "PercentFormatExpectedSequence": "F503",
    "PercentFormatExtraNamedArguments": "F504",
    "PercentFormatMissingArgument": "F505",
    "PercentFormatMixedPositionalAndNamed": "F506",
    "PercentFormatPositionalCountMismatch": "F507",
    "PercentFormatStarRequiresSequence": "F508",
    "PercentFormatUnsupportedFormatCharacter": "F509",
    "StringDotFormatInvalidFormat": "F521",
    "StringDotFormatExtraNamedArguments": "F522",
    "StringDotFormatExtraPositionalArguments": "F523",
    "StringDotFormatMissingArgument": "F524",
    "StringDotFormatMixingAutomatic": "F525",
    "FStringMissingPlaceholders": "F541",
    "MultiValueRepeatedKeyLiteral": "F601",
    "MultiValueRepeatedKeyVariable": "F602",
    "TooManyExpressionsInStarredAssignment": "F621",
    "TwoStarredExpressions": "F622",
    "AssertTuple": "F631",
    "IsLiteral": "F632",
    "InvalidPrintSyntax": "F633",
    "IfTuple": "F634",
    "BreakOutsideLoop": "F701",
    "ContinueOutsideLoop": "F702",
    "YieldOutsideFunction": "F704",
    "ReturnOutsideFunction": "F706",
    "DefaultExceptNotLast": "F707",
    "DoctestSyntaxError": "F721",
    "ForwardAnnotationSyntaxError": "F722",
    "RedefinedWhileUnused": "F811",
    "UndefinedName": "F821",
    "UndefinedExport": "F822",
    "UndefinedLocal": "F823",
    "DuplicateArgument": "F831",
    "UnusedVariable": "F841",
    "UnusedAnnotation": "F842",
    "RaiseNotImplemented": "F901",
}


class FlakesChecker(pyflakes.checker.Checker):
    """Subclass the Pyflakes checker to conform with the flake8 API."""

    with_doctest = False
    include_in_doctest: list[str] = []
    exclude_from_doctest: list[str] = []

    def __init__(self, tree: ast.AST, filename: str) -> None:
        """Initialize the PyFlakes plugin with an AST tree and filename."""
        filename = utils.normalize_path(filename)
        with_doctest = self.with_doctest
        included_by = [
            include
            for include in self.include_in_doctest
            if include != "" and filename.startswith(include)
        ]
        if included_by:
            with_doctest = True

        for exclude in self.exclude_from_doctest:
            if exclude != "" and filename.startswith(exclude):
                with_doctest = False
                overlapped_by = [
                    include
                    for include in included_by
                    if include.startswith(exclude)
                ]

                if overlapped_by:
                    with_doctest = True

        super().__init__(tree, filename=filename, withDoctest=with_doctest)

    @classmethod
    def add_options(cls, parser: OptionManager) -> None:
        """Register options for PyFlakes on the Flake8 OptionManager."""
        parser.add_option(
            "--builtins",
            parse_from_config=True,
            comma_separated_list=True,
            help="define more built-ins, comma separated",
        )
        parser.add_option(
            "--doctests",
            default=False,
            action="store_true",
            parse_from_config=True,
            help="also check syntax of the doctests",
        )
        parser.add_option(
            "--include-in-doctest",
            default="",
            dest="include_in_doctest",
            parse_from_config=True,
            comma_separated_list=True,
            normalize_paths=True,
            help="Run doctests only on these files",
        )
        parser.add_option(
            "--exclude-from-doctest",
            default="",
            dest="exclude_from_doctest",
            parse_from_config=True,
            comma_separated_list=True,
            normalize_paths=True,
            help="Skip these files when running doctests",
        )

    @classmethod
    def parse_options(cls, options: argparse.Namespace) -> None:
        """Parse option values from Flake8's OptionManager."""
        if options.builtins:
            cls.builtIns = cls.builtIns.union(options.builtins)
        cls.with_doctest = options.doctests

        if options.include_in_doctest or options.exclude_from_doctest:
            LOG.warning(
                "--include-in-doctest / --exclude-from-doctest will be "
                "removed in a future version.  see PyCQA/flake8#1747"
            )

        included_files = []
        for included_file in options.include_in_doctest:
            if included_file == "":
                continue
            if not included_file.startswith((os.sep, "./", "~/")):
                included_files.append(f"./{included_file}")
            else:
                included_files.append(included_file)
        cls.include_in_doctest = utils.normalize_paths(included_files)

        excluded_files = []
        for excluded_file in options.exclude_from_doctest:
            if excluded_file == "":
                continue
            if not excluded_file.startswith((os.sep, "./", "~/")):
                excluded_files.append(f"./{excluded_file}")
            else:
                excluded_files.append(excluded_file)
        cls.exclude_from_doctest = utils.normalize_paths(excluded_files)

        inc_exc = set(cls.include_in_doctest).intersection(
            cls.exclude_from_doctest
        )
        if inc_exc:
            raise ValueError(
                f"{inc_exc!r} was specified in both the "
                f"include-in-doctest and exclude-from-doctest "
                f"options. You are not allowed to specify it in "
                f"both for doctesting."
            )

    def run(self) -> Generator[tuple[int, int, str, type[Any]], None, None]:
        """Run the plugin."""
        for message in self.messages:
            col = getattr(message, "col", 0)
            yield (
                message.lineno,
                col,
                "{} {}".format(
                    FLAKE8_PYFLAKES_CODES.get(type(message).__name__, "F999"),
                    message.message % message.message_args,
                ),
                message.__class__,
            )
