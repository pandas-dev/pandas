import platform
import unittest
from unittest import skipUnless
from unittest.mock import NonCallableMock
from itertools import chain
from datetime import datetime
from contextlib import redirect_stdout
from io import StringIO

from numba.tests.support import TestCase
import numba.misc.numba_sysinfo as nsi


class TestSysInfo(TestCase):

    def setUp(self):
        super(TestSysInfo, self).setUp()
        self.info = nsi.get_sysinfo()
        self.safe_contents = {
            int: (
                nsi._cpu_count,
            ),
            float: (
                nsi._runtime,
            ),
            str: (
                nsi._machine,
                nsi._cpu_name,
                nsi._platform_name,
                nsi._os_name,
                nsi._os_version,
                nsi._python_comp,
                nsi._python_impl,
                nsi._python_version,
                nsi._llvm_version,
                nsi._numpy_version,
            ),
            bool: (
                nsi._cu_dev_init,
                nsi._svml_state,
                nsi._svml_loaded,
                nsi._svml_operational,
                nsi._llvm_svml_patched,
                nsi._tbb_thread,
                nsi._openmp_thread,
                nsi._wkq_thread,
                nsi._numpy_AVX512_SKX_detected,
            ),
            list: (
                nsi._errors,
                nsi._warnings,
            ),
            dict: (
                nsi._numba_env_vars,
            ),
            datetime: (
                nsi._start,
                nsi._start_utc,
            ),
        }
        self.safe_keys = chain(*self.safe_contents.values())

    def tearDown(self):
        super(TestSysInfo, self).tearDown()
        # System info might contain long strings or lists so delete it.
        del self.info

    def test_has_safe_keys(self):
        for k in self.safe_keys:
            with self.subTest(k=k):
                self.assertIn(k, self.info)

    def test_safe_content_type(self):
        for t, keys in self.safe_contents.items():
            for k in keys:
                with self.subTest(k=k):
                    self.assertIsInstance(self.info[k], t)

    def test_has_no_error(self):
        self.assertFalse(self.info[nsi._errors])

    def test_display_empty_info(self):
        output = StringIO()
        with redirect_stdout(output):
            res = nsi.display_sysinfo({})
        self.assertIsNone(res)
        output.close()


class TestSysInfoWithPsutil(TestCase):

    mem_total = 2 * 1024 ** 2  # 2_097_152
    mem_available = 1024 ** 2  # 1_048_576
    cpus_list = [1, 2]

    def setUp(self):
        super(TestSysInfoWithPsutil, self).setUp()
        self.psutil_orig_state = nsi._psutil_import
        # Mocking psutil
        nsi._psutil_import = True
        nsi.psutil = NonCallableMock()
        vm = nsi.psutil.virtual_memory.return_value
        vm.total = self.mem_total
        vm.available = self.mem_available
        if platform.system() in ('Linux', 'Windows',):
            # cpu_affiniy only available on Linux and Windows
            proc = nsi.psutil.Process.return_value
            proc.cpu_affinity.return_value = self.cpus_list
        else:
            nsi.psutil.Process.return_value = None

        self.info = nsi.get_os_spec_info(platform.system())

    def tearDown(self):
        super(TestSysInfoWithPsutil, self).tearDown()
        nsi._psutil_import = self.psutil_orig_state

    def test_has_all_data(self):
        keys = (nsi._mem_total, nsi._mem_available)
        for k in keys:
            with self.subTest(k=k):
                self.assertIn(k, self.info.keys())
                self.assertIsInstance(self.info[k], int)

    def test_has_correct_values(self):
        self.assertEqual(self.info[nsi._mem_total], self.mem_total)
        self.assertEqual(self.info[nsi._mem_available], self.mem_available)

    @skipUnless(platform.system() in ('Linux', 'Windows'),
                "CPUs allowed info only available on Linux and Windows")
    def test_cpus_list(self):
        self.assertEqual(self.info[nsi._cpus_allowed], len(self.cpus_list))
        self.assertEqual(self.info[nsi._cpus_list],
                         ' '.join(str(n) for n in self.cpus_list))


class TestSysInfoWithoutPsutil(TestCase):

    def setUp(self):
        super(TestSysInfoWithoutPsutil, self).setUp()
        self.psutil_orig_state = nsi._psutil_import
        nsi._psutil_import = False
        self.info = nsi.get_os_spec_info(platform.system())

    def tearDown(self):
        super(TestSysInfoWithoutPsutil, self).tearDown()
        nsi._psutil_import = self.psutil_orig_state

    def test_has_all_data(self):
        keys = (nsi._mem_total, nsi._mem_available)
        for k in keys:
            with self.subTest(k=k):
                self.assertIn(k, self.info.keys())
                self.assertIsInstance(self.info[k], int)


class TestPlatformSpecificInfo(TestCase):

    def setUp(self):
        self.plat_spec_info = {
            'Linux': {
                str: (nsi._libc_version,),
            },
            'Windows': {
                str: (nsi._os_spec_version,),
            },
            'Darwin': {
                str: (nsi._os_spec_version,),
            },
        }
        self.os_name = platform.system()
        self.contents = self.plat_spec_info.get(self.os_name, {})
        self.info = nsi.get_os_spec_info(self.os_name)

    def test_has_all_data(self):
        keys = chain(*self.contents.values())
        for k in keys:
            with self.subTest(k=k):
                self.assertIn(k, self.info.keys())

    def test_content_type(self):
        for t, keys in self.contents.items():
            for k in keys:
                with self.subTest(k=k):
                    self.assertIsInstance(self.info[k], t)


if __name__ == '__main__':
    unittest.main()
