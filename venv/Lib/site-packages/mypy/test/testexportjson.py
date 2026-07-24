"""Test cases for the mypy cache JSON export tool."""

from __future__ import annotations

import json
import os
import re
import sys

from mypy import build
from mypy.errors import CompileError
from mypy.exportjson import convert_binary_cache_meta_to_json, convert_binary_cache_to_json
from mypy.modulefinder import BuildSource
from mypy.options import Options
from mypy.test.config import test_temp_dir
from mypy.test.data import DataDrivenTestCase, DataSuite
from mypy.test.helpers import assert_string_arrays_equal


class TypeExportSuite(DataSuite):
    required_out_section = True
    files = ["exportjson.test"]

    def run_case(self, testcase: DataDrivenTestCase) -> None:
        error = False
        src = "\n".join(testcase.input)
        is_meta = testcase.name.endswith("_meta")
        try:
            options = Options()
            options.use_builtins_fixtures = True
            options.show_traceback = True
            options.allow_empty_bodies = True
            options.fixed_format_cache = True
            options.sqlite_cache = False
            fnam = os.path.join(self.base_path, "main.py")
            with open(fnam, "w") as f:
                f.write(src)
            result = build.build(
                sources=[BuildSource(fnam, "main")], options=options, alt_lib_path=test_temp_dir
            )
            a = result.errors
            error = bool(a)

            major, minor = sys.version_info[:2]
            cache_dir = os.path.join(".mypy_cache", f"{major}.{minor}")

            for module in result.files:
                if module in (
                    "builtins",
                    "typing",
                    "_typeshed",
                    "__future__",
                    "typing_extensions",
                    "sys",
                    "collections",
                ):
                    continue
                fnam = os.path.join(cache_dir, f"{module}.data.ff")
                if not is_meta:
                    with open(fnam, "rb") as f:
                        json_data = convert_binary_cache_to_json(f.read(), implicit_names=False)
                else:
                    meta_fnam = os.path.join(cache_dir, f"{module}.meta.ff")
                    with open(meta_fnam, "rb") as f:
                        json_data = convert_binary_cache_meta_to_json(f.read(), fnam)
                for line in json.dumps(json_data, indent=4).splitlines():
                    if '"path": ' in line:
                        # The source file path is unpredictable, so filter it out
                        line = re.sub(r'"[^"]+\.pyi?"', "...", line)
                    if is_meta:
                        if '"version_id"' in line:
                            line = re.sub(r'"[0-9][^"]+"', "...", line)
                        if '"mtime"' in line or '"data_mtime"' in line:
                            line = re.sub(r": [0-9]+", ": ...", line)
                        if '"platform"' in line:
                            line = re.sub(': "[^"]+"', ": ...", line)
                        if '"hash"' not in line:
                            # Some hashes are unpredictable so filter them out
                            line = re.sub(r'"[a-f0-9]{40}"', '"<hash>"', line)
                    assert "ERROR" not in line, line
                    a.append(line)
        except CompileError as e:
            a = e.messages
            error = True
        if error or "\n".join(testcase.output).strip() != "<not checked>":
            if is_meta and sys.platform == "win32":
                out = filter_platform_specific(testcase.output)
                a = filter_platform_specific(a)
            else:
                out = testcase.output
            assert_string_arrays_equal(
                out, a, f"Invalid output ({testcase.file}, line {testcase.line})"
            )


def filter_platform_specific(lines: list[str]) -> list[str]:
    return [l for l in lines if '"size":' not in l and '"hash":' not in l]
