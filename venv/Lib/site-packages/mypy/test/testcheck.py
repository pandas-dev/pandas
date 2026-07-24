"""Type checker test cases"""

from __future__ import annotations

import os
import re
import sys
import sysconfig
import tempfile
from pathlib import Path

from mypy import build
from mypy.errors import CompileError
from mypy.modulefinder import BuildSource, FindModuleCache, SearchPaths
from mypy.test.config import test_data_prefix, test_temp_dir
from mypy.test.data import DataDrivenTestCase, DataSuite, FileOperation, module_from_path
from mypy.test.helpers import (
    assert_module_equivalence,
    assert_string_arrays_equal,
    assert_target_equivalence,
    check_test_output_files,
    find_test_files,
    normalize_error_messages,
    parse_options,
    perform_file_operations,
    remove_typevar_ids,
)
from mypy.test.update_data import update_testcase_output

try:
    if sys.version_info >= (3, 14) and bool(sysconfig.get_config_var("Py_GIL_DISABLED")):
        # lxml doesn't support free-threading yet
        lxml = None
    else:
        import lxml  # type: ignore[import-untyped]
except ImportError:
    lxml = None


import pytest

# List of files that contain test case descriptions.
# Includes all check-* files with the .test extension in the test-data/unit directory
typecheck_files = find_test_files(pattern="check-*.test")

# Tests that use Python version specific features:
if sys.version_info < (3, 11):
    typecheck_files.remove("check-python311.test")
if sys.version_info < (3, 12):
    typecheck_files.remove("check-python312.test")
if sys.version_info < (3, 13):
    typecheck_files.remove("check-python313.test")
if sys.version_info < (3, 14):
    typecheck_files.remove("check-python314.test")


class TypeCheckSuite(DataSuite):
    files = typecheck_files

    def run_case(self, testcase: DataDrivenTestCase) -> None:
        if os.path.basename(testcase.file) == "check-modules-case.test":
            with tempfile.NamedTemporaryFile(prefix="test", dir=".") as temp_file:
                temp_path = Path(temp_file.name)
                if not temp_path.with_name(temp_path.name.upper()).exists():
                    pytest.skip("File system is not case-insensitive")
        if lxml is None and os.path.basename(testcase.file) == "check-reports.test":
            pytest.skip("Cannot import lxml. Is it installed?")
        incremental = (
            "incremental" in testcase.name.lower()
            or "incremental" in testcase.file
            or "serialize" in testcase.file
        )
        if incremental:
            # Incremental tests are run once with a cold cache, once with a warm cache.
            # Expect success on first run, errors from testcase.output (if any) on second run.
            num_steps = max([2] + list(testcase.output2.keys()))
            # Check that there are no file changes beyond the last run (they would be ignored).
            for dn, dirs, files in os.walk(os.curdir):
                for file in files:
                    m = re.search(r"\.([2-9])$", file)
                    if m and int(m.group(1)) > num_steps:
                        raise ValueError(
                            "Output file {} exists though test case only has {} runs".format(
                                file, num_steps
                            )
                        )
            steps = testcase.find_steps()
            for step in range(1, num_steps + 1):
                idx = step - 2
                ops = steps[idx] if idx < len(steps) and idx >= 0 else []
                self.run_case_once(testcase, ops, step, True)
        else:
            self.run_case_once(testcase)

    def _sort_output_if_needed(self, testcase: DataDrivenTestCase, a: list[str]) -> None:
        idx = testcase.output_inline_start
        if not testcase.files or idx == len(testcase.output):
            return

        def _filename(_msg: str) -> str:
            return _msg.partition(":")[0]

        file_weights = {file: idx for idx, file in enumerate(_filename(msg) for msg in a)}
        testcase.output[idx:] = sorted(
            testcase.output[idx:], key=lambda msg: file_weights.get(_filename(msg), -1)
        )

    def run_case_once(
        self,
        testcase: DataDrivenTestCase,
        operations: list[FileOperation] | None = None,
        incremental_step: int = 0,
        is_incremental: bool = False,
    ) -> None:
        if operations is None:
            operations = []
        original_program_text = "\n".join(testcase.input)
        module_data = self.parse_module(original_program_text, incremental_step)

        # Unload already loaded plugins, they may be updated.
        for file, _ in testcase.files:
            module = module_from_path(file)
            if module.endswith("_plugin") and module in sys.modules:
                del sys.modules[module]
        if incremental_step == 0 or incremental_step == 1:
            # In run 1, copy program text to program file.
            for module_name, program_path, program_text in module_data:
                if module_name == "__main__":
                    with open(program_path, "w", encoding="utf8") as f:
                        f.write(program_text)
                    break
        elif incremental_step > 1:
            # In runs 2+, copy *.[num] files to * files.
            perform_file_operations(operations)

        # Parse options after moving files (in case mypy.ini is being moved).
        options = parse_options(original_program_text, testcase, incremental_step)
        options.use_builtins_fixtures = True
        options.show_traceback = True
        options.native_parser = bool(os.environ.get("TEST_NATIVE_PARSER"))
        options.reveal_verbose_types = not testcase.name.endswith("_no_verbose_reveal")

        if options.num_workers:
            options.fixed_format_cache = True
            options.native_parser = True
            if testcase.output_files:
                raise pytest.skip("Reports are not supported in parallel mode")
            # Note: do not use this unless really needed!
            if testcase.name.endswith("_no_parallel"):
                raise pytest.skip("Test not supported in parallel mode yet")
        else:
            if testcase.name.endswith("_parallel_only"):
                raise pytest.skip("Test is only for parallel mode")

        if options.native_parser and testcase.name.endswith("_no_native_parse"):
            raise pytest.skip("Test not supported by native parser yet")

        # Enable some options automatically based on test file name.
        if "columns" in testcase.file:
            options.show_column_numbers = True
        if "errorcodes" in testcase.file:
            options.hide_error_codes = False
        if "abstract" not in testcase.file:
            options.allow_empty_bodies = not testcase.name.endswith("_no_empty")

        if incremental_step and options.incremental or options.num_workers > 0:
            # Don't overwrite # flags: --no-incremental in incremental test cases
            options.incremental = True
        else:
            options.incremental = False
            # Don't waste time writing cache unless we are specifically looking for it
            if not testcase.writescache:
                options.cache_dir = os.devnull

        sources = []
        for module_name, program_path, program_text in module_data:
            # Always set to None, so we're forced to reread the module in incremental mode
            sources.append(
                BuildSource(program_path, module_name, None if incremental_step else program_text)
            )

        plugin_dir = os.path.join(test_data_prefix, "plugins")

        worker_env = None
        if options.num_workers > 0:
            worker_env = os.environ.copy()
            # Make sure we are running tests with current worktree files, *not* with
            # an installed version of mypy.
            root_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
            worker_env["PYTHONPATH"] = os.pathsep.join([root_dir, plugin_dir])
            worker_env["MYPY_TEST_PREFIX"] = root_dir
            worker_env["MYPY_ALT_LIB_PATH"] = test_temp_dir

        sys.path.insert(0, plugin_dir)
        res = None
        blocker = False
        try:
            res = build.build(
                sources=sources, options=options, alt_lib_path=test_temp_dir, worker_env=worker_env
            )
            a = res.errors
        except CompileError as e:
            a = e.messages
            blocker = True
        finally:
            assert sys.path[0] == plugin_dir
            del sys.path[0]

        if testcase.normalize_output:
            a = normalize_error_messages(a)

        # Make sure error messages match
        if incremental_step < 2:
            if incremental_step == 1:
                msg = "Unexpected type checker output in incremental, run 1 ({}, line {})"
            else:
                assert incremental_step == 0
                msg = "Unexpected type checker output ({}, line {})"
            self._sort_output_if_needed(testcase, a)
            output = testcase.output
        else:
            msg = (
                f"Unexpected type checker output in incremental, run {incremental_step}"
                + " ({}, line {})"
            )
            output = testcase.output2.get(incremental_step, [])

        if output != a and testcase.config.getoption("--update-data", False):
            update_testcase_output(testcase, a, incremental_step=incremental_step)

        ignore_modules_order = False
        if options.num_workers > 0:
            ignore_modules_order = is_incremental
            # TypeVarIds are not stable in parallel checking, normalize.
            a = remove_typevar_ids(a)
            output = remove_typevar_ids(output)
        assert_string_arrays_equal(
            output,
            a,
            msg.format(testcase.file, testcase.line),
            ignore_modules_order=ignore_modules_order,
        )

        if res:
            if options.cache_dir != os.devnull:
                self.verify_cache(module_data, res.manager, blocker, incremental_step)

            name = "targets"
            if incremental_step:
                name += str(incremental_step + 1)
            expected = testcase.expected_fine_grained_targets.get(incremental_step + 1)
            actual = [
                target
                for module, target in res.manager.processed_targets
                if module in testcase.test_modules
            ]
            # TODO: check targets in parallel mode (e.g. per SCC).
            if options.num_workers == 0 and expected is not None:
                assert_target_equivalence(name, expected, actual)
            if incremental_step > 1:
                suffix = "" if incremental_step == 2 else str(incremental_step - 1)
                expected_rechecked = testcase.expected_rechecked_modules.get(incremental_step - 1)
                if expected_rechecked is not None:
                    assert_module_equivalence(
                        "rechecked" + suffix, expected_rechecked, res.manager.rechecked_modules
                    )
                expected_stale = testcase.expected_stale_modules.get(incremental_step - 1)
                if expected_stale is not None:
                    assert_module_equivalence(
                        "stale" + suffix, expected_stale, res.manager.stale_modules
                    )
            res.manager.metastore.close()

        if testcase.output_files:
            check_test_output_files(testcase, incremental_step, strip_prefix="tmp/")

    def verify_cache(
        self,
        module_data: list[tuple[str, str, str]],
        manager: build.BuildManager,
        blocker: bool,
        step: int,
    ) -> None:
        if not blocker:
            # There should be valid cache metadata for each module except
            # in case of a blocking error in themselves or one of their
            # dependencies.
            modules = self.find_module_files(manager)
            modules.update({module_name: path for module_name, path, text in module_data})
            missing_paths = self.find_missing_cache_files(modules, manager)
            if missing_paths:
                raise AssertionError(f"cache data missing for {missing_paths} on run {step}")
        assert os.path.isfile(os.path.join(manager.options.cache_dir, ".gitignore"))
        cachedir_tag = os.path.join(manager.options.cache_dir, "CACHEDIR.TAG")
        assert os.path.isfile(cachedir_tag)
        with open(cachedir_tag) as f:
            assert f.read().startswith("Signature: 8a477f597d28d172789f06886806bc55")

    def find_module_files(self, manager: build.BuildManager) -> dict[str, str]:
        return {id: module.path for id, module in manager.modules.items()}

    def find_missing_cache_files(
        self, modules: dict[str, str], manager: build.BuildManager
    ) -> set[str]:
        ignore_errors = True
        missing = {}
        for id, path in modules.items():
            meta_pair = build.find_cache_meta(id, path, manager)
            if meta_pair is None:
                missing[id] = path
            elif not build.validate_meta(meta_pair[0], id, path, ignore_errors, manager):
                missing[id] = path
        return set(missing.values())

    def parse_module(
        self, program_text: str, incremental_step: int = 0
    ) -> list[tuple[str, str, str]]:
        """Return the module and program names for a test case.

        Normally, the unit tests will parse the default ('__main__')
        module and follow all the imports listed there. You can override
        this behavior and instruct the tests to check multiple modules
        by using a comment like this in the test case input:

          # cmd: mypy -m foo.bar foo.baz

        You can also use `# cmdN:` to have a different cmd for incremental
        step N (2, 3, ...).

        Return a list of tuples (module name, file name, program text).
        """
        m = re.search("# cmd: mypy -m ([a-zA-Z0-9_. ]+)$", program_text, flags=re.MULTILINE)
        if incremental_step > 1:
            alt_regex = f"# cmd{incremental_step}: mypy -m ([a-zA-Z0-9_. ]+)$"
            alt_m = re.search(alt_regex, program_text, flags=re.MULTILINE)
            if alt_m is not None:
                # Optionally return a different command if in a later step
                # of incremental mode, otherwise default to reusing the
                # original cmd.
                m = alt_m

        if m:
            # The test case wants to use a non-default main
            # module. Look up the module and give it as the thing to
            # analyze.
            module_names = m.group(1)
            out = []
            search_paths = SearchPaths((test_temp_dir,), (), (), ())
            cache = FindModuleCache(search_paths, fscache=None, options=None)
            for module_name in module_names.split(" "):
                path = cache.find_module(module_name)
                assert isinstance(path, str), f"Can't find ad hoc case file: {module_name}"
                with open(path, encoding="utf8") as f:
                    program_text = f.read()
                out.append((module_name, path, program_text))
            return out
        else:
            return [("__main__", "main", program_text)]
