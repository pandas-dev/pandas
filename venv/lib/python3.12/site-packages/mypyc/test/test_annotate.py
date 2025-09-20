"""Test cases for annotating source code to highlight inefficiencies."""

from __future__ import annotations

import os.path

from mypy.errors import CompileError
from mypy.test.config import test_temp_dir
from mypy.test.data import DataDrivenTestCase
from mypyc.annotate import generate_annotations, get_max_prio
from mypyc.ir.pprint import format_func
from mypyc.test.testutil import (
    ICODE_GEN_BUILTINS,
    MypycDataSuite,
    assert_test_output,
    build_ir_for_single_file2,
    infer_ir_build_options_from_test_name,
    remove_comment_lines,
    use_custom_builtins,
)

files = ["annotate-basic.test"]


class TestReport(MypycDataSuite):
    files = files
    base_path = test_temp_dir
    optional_out = True

    def run_case(self, testcase: DataDrivenTestCase) -> None:
        """Perform a runtime checking transformation test case."""
        options = infer_ir_build_options_from_test_name(testcase.name)
        if options is None:
            # Skipped test case
            return
        with use_custom_builtins(os.path.join(self.data_prefix, ICODE_GEN_BUILTINS), testcase):
            expected_output = remove_comment_lines(testcase.output)

            # Parse "# A: <message>" comments.
            for i, line in enumerate(testcase.input):
                if "# A:" in line:
                    msg = line.rpartition("# A:")[2].strip()
                    expected_output.append(f"main:{i + 1}: {msg}")

            ir = None
            try:
                ir, tree, type_map, mapper = build_ir_for_single_file2(testcase.input, options)
            except CompileError as e:
                actual = e.messages
            else:
                annotations = generate_annotations("native.py", tree, ir, type_map, mapper)
                actual = []
                for line_num, line_anns in sorted(
                    annotations.annotations.items(), key=lambda it: it[0]
                ):
                    anns = get_max_prio(line_anns)
                    str_anns = [a.message for a in anns]
                    s = " ".join(str_anns)
                    actual.append(f"main:{line_num}: {s}")

            try:
                assert_test_output(testcase, actual, "Invalid source code output", expected_output)
            except BaseException:
                if ir:
                    print("Generated IR:\n")
                    for fn in ir.functions:
                        if fn.name == "__top_level__":
                            continue
                        for s in format_func(fn):
                            print(s)
                raise
