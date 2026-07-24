from __future__ import annotations

import getopt
import os
import platform
import re
import shutil
import sys
import traceback
import unittest
from collections.abc import Callable
from textwrap import dedent
from types import ModuleType
from typing import NoReturn

import win32com

this_folder = os.path.dirname(__file__)
win32com_src_dir = os.path.abspath(os.path.join(this_folder, ".."))
# We'd prefer the win32com namespace to be the parent of __file__ - ie, our source-tree,
# rather than the version installed - otherwise every .py change needs a full install to
# test!
# We can't patch win32comext as most of them have a .pyd in their root :(
# This clearly isn't ideal or perfect :)
win32com.__path__[0] = win32com_src_dir

import pythoncom
import win32com.client.gencache
from pywin32_testutil import TestLoader, TestRunner
from win32com.test.util import (
    CapturingFunctionTestCase,
    CheckClean,
    RegisterPythonServer,
    ShellTestCase,
    TestCase,
)

verbosity = 1  # default unittest verbosity.


def CleanGenerated() -> None:
    if os.path.isdir(win32com.__gen_path__):
        if verbosity > 1:
            print(f"Deleting files from", win32com.__gen_path__)
        shutil.rmtree(win32com.__gen_path__)

    win32com.client.gencache.__init__()  # Reset


def RemoveRefCountOutput(data: str) -> str:
    while True:
        last_line_pos = data.rfind("\n")
        if not re.match(r"\[\d+ refs\]", data[last_line_pos + 1 :]):
            break
        if last_line_pos < 0:
            # All the output
            return ""
        data = data[:last_line_pos]

    return data


def ExecuteSilentlyIfOK(cmd: str, testcase: TestCase) -> str:
    f = os.popen(cmd)
    data = f.read().strip()
    rc = f.close()
    if rc:
        print(data)
        testcase.fail(f"Executing '{cmd}' failed ({rc})")
    # for "_d" builds, strip the '[xxx refs]' line
    return RemoveRefCountOutput(data)


@unittest.skipIf(
    platform.machine() == "ARM64",
    "PyCOMTest cannot currently be run on ARM64 "
    + "due to lacking win32com.universal implementation "
    + "in com/win32com/src/univgw.cpp",
)
class PyCOMTest(TestCase):
    no_leak_tests = True  # done by the test itself

    def testit(self) -> None:
        # Check that the item is registered, so we get the correct
        # 'skipped' behaviour (and recorded as such) rather than either
        # error or silence due to non-registration.
        RegisterPythonServer(
            os.path.join(win32com_src_dir, "servers", "test_pycomtest.py"),
            "Python.Test.PyCOMTest",
        )

        # Execute testPyComTest in its own process so it can play
        # with the Python thread state
        fname = os.path.join(this_folder, "testPyComTest.py")
        cmd = f'{sys.executable} "{fname}" -q 2>&1'
        data = ExecuteSilentlyIfOK(cmd, self)


class PippoTest(TestCase):
    def testit(self) -> None:
        # Check we are registered before spawning the process.
        from win32com.test import pippo_server

        RegisterPythonServer(pippo_server.__file__, "Python.Test.Pippo")

        python = sys.executable
        fname = os.path.join(this_folder, "testPippo.py")
        cmd = f'{python} "{fname}" 2>&1'
        ExecuteSilentlyIfOK(cmd, self)


# This is a list of "win32com.test.???" module names, optionally with a
# function in that module if the module isn't unitest based...
unittest_modules = [
    # Level 1 tests - fast and few dependencies - good for CI!
    [
        "testIterators",
        "testvbscript_regexp",
        "testStorage",
        "testStreams",
        "testWMI",
        "policySemantics",
        "testShell",
        "testROT",
        "testxslt",
        "testCollections",
        "errorSemantics.test",
        "testArrays",
        "testClipboard",
        "testConversionErrors",
    ],
    # Level 2 tests - wants our demo COM objects registered.
    # (these are strange; on GitHub CI they get further than expected when
    # our objects are not installed, so fail to quietly fail with "can't
    # register" like they do locally. So really just a nod to CI)
    ["testAXScript", "testDictionary", "testServers", "testvb", "testMarshal"],
    # Level 3 tests - Requires Office or other non-free stuff.
    [
        "testMSOffice.TestAll",
        "testMSOfficeEvents.test",
        "testAccess.test",
        "testExplorer.TestAll",
        "testExchange.test",
    ],
    # Level 4 tests - we try and run `makepy` over every typelib installed!
    ["testmakepy.TestAll"],
]

# A list of other unittest modules we use - these are fully qualified module
# names and the module is assumed to be unittest based.
unittest_other_modules = [
    # Level 1 tests.
    ["win32com.directsound.test.ds_test"],
    # Level 2 tests.
    [],
    # Level 3 tests.
    [],
    # Level 4 tests.
    [],
]


output_checked_programs = [
    # Level 1 tests.
    [],
    # Level 2 tests.
    [
        ("cscript.exe /nologo //E:vbscript testInterp.vbs", "VBScript test worked OK"),
        (
            "cscript.exe /nologo //E:vbscript testDictionary.vbs",
            "VBScript has successfully tested Python.Dictionary",
        ),
    ],
    # Level 3 tests
    [],
    # Level 4 tests.
    [],
]

custom_test_cases = [
    # Level 1 tests.
    [
        PyCOMTest,
    ],
    # Level 2 tests.
    [
        PippoTest,
    ],
    # Level 3 tests
    [],
    # Level 4 tests.
    [],
]


def get_test_mod_and_func(
    test_name: str, import_failures: list[tuple[str, BaseException]]
) -> tuple[None, None] | tuple[ModuleType, Callable[[], object] | None]:
    if test_name.find(".") > 0:
        mod_name, func_name = test_name.split(".")
    else:
        mod_name = test_name
        func_name = None
    fq_mod_name = "win32com.test." + mod_name
    try:
        __import__(fq_mod_name)
        mod = sys.modules[fq_mod_name]
    except Exception as error:
        import_failures.append((mod_name, error))
        return None, None
    func = None if func_name is None else getattr(mod, func_name)
    return mod, func


# Return a test suite all loaded with the tests we want to run
def make_test_suite(
    test_level: int,
) -> tuple[unittest.TestSuite, list[tuple[str, BaseException]]]:
    suite = unittest.TestSuite()
    import_failures: list[tuple[str, BaseException]] = []
    loader = TestLoader()
    for i in range(test_level):
        for mod_name in unittest_modules[i]:
            mod, func = get_test_mod_and_func(mod_name, import_failures)
            if mod is None:
                raise ModuleNotFoundError(f"no such module '{mod_name}'")
            if func is not None:
                test = CapturingFunctionTestCase(func, description=mod_name)
            else:
                if hasattr(mod, "suite"):
                    test = mod.suite()
                else:
                    test = loader.loadTestsFromModule(mod)
            assert test.countTestCases() > 0, f"No tests loaded from {mod!r}"
            suite.addTest(test)
        for cmd, output in output_checked_programs[i]:
            suite.addTest(ShellTestCase(cmd, output))

        for test_class in custom_test_cases[i]:
            suite.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(test_class))
    # other "normal" unittest modules.
    for i in range(test_level):
        for mod_name in unittest_other_modules[i]:
            try:
                __import__(mod_name)
            except Exception as error:
                import_failures.append((mod_name, error))
                continue

            mod = sys.modules[mod_name]
            if hasattr(mod, "suite"):
                test = mod.suite()
            else:
                test = loader.loadTestsFromModule(mod)
            assert test.countTestCases() > 0, f"No tests loaded from {mod!r}"
            suite.addTest(test)

    return suite, import_failures


def usage(why: str) -> NoReturn:
    print(f"""\
{why}

win32com test suite
usage: testall [-v] test_level
  where test_level is an integer 1-4.  Level 1 tests are quick,
  level 2 tests invoke Word, IE etc, level 3 take ages!
  level 4 we try and run `makepy` over every typelib installed""")
    sys.exit(1)


if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "v")
    except getopt.error as why:
        usage(why)
    for opt, val in opts:
        if opt == "-v":
            verbosity += 1

    if len(args) > 1:
        usage(
            "Only a single test level argument is supported. "
            + "Test names are not supported yet"
        )
    try:
        test_level = int(args[0])
    except IndexError:
        test_level = 2  # default to quick test with local objects
    except ValueError:
        usage("Test names are not supported yet")
    if test_level < 1 or test_level > 4:
        usage("Only levels 1-4 are supported")
    CleanGenerated()

    suite, import_failures = make_test_suite(test_level)
    if verbosity:
        if hasattr(sys, "gettotalrefcount"):
            print(
                dedent("""\
                This is a debug build - memory leak tests will also be run.
                These tests may take *many* minutes to run - be patient!
                (running from python.exe will avoid these leak tests""")
            )
        print(
            f"Executing level {test_level} tests - {suite.countTestCases()} test cases will be run"
        )
        if verbosity == 1 and suite.countTestCases() < 70:
            # A little row of markers so the dots show how close to finished
            print("|" * suite.countTestCases())
    testResult = TestRunner(verbosity=verbosity).run(suite)
    if import_failures:
        testResult.stream.writeln(
            "*** The following test modules could not be imported ***"
        )
        for mod_name, error in import_failures:
            desc = "\n".join(traceback.format_exception_only(type(error), error))
            testResult.stream.write(f"{mod_name}: {desc}")
        testResult.stream.writeln(
            f"*** {len(import_failures)} test(s) could not be run ***"
        )

    # re-print unit-test error here so it is noticed
    if not testResult.wasSuccessful():
        print("*" * 20, "- unittest tests FAILED")

    CheckClean()
    pythoncom.CoUninitialize()
    CleanGenerated()
    if not testResult.wasSuccessful():
        sys.exit(1)
