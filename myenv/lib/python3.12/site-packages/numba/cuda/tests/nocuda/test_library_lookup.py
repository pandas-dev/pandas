import sys
import os
import multiprocessing as mp
import warnings

from numba.core.config import IS_WIN32, IS_OSX
from numba.core.errors import NumbaWarning
from numba.cuda.cudadrv import nvvm
from numba.cuda.testing import (
    unittest,
    skip_on_cudasim,
    SerialMixin,
    skip_unless_conda_cudatoolkit,
)
from numba.cuda.cuda_paths import (
    _get_libdevice_path_decision,
    _get_nvvm_path_decision,
    _get_cudalib_dir_path_decision,
    get_system_ctk,
)


has_cuda = nvvm.is_available()
has_mp_get_context = hasattr(mp, 'get_context')


class LibraryLookupBase(SerialMixin, unittest.TestCase):
    def setUp(self):
        ctx = mp.get_context('spawn')

        qrecv = ctx.Queue()
        qsend = ctx.Queue()
        self.qsend = qsend
        self.qrecv = qrecv
        self.child_process = ctx.Process(
            target=check_lib_lookup,
            args=(qrecv, qsend),
            daemon=True,
        )
        self.child_process.start()

    def tearDown(self):
        self.qsend.put(self.do_terminate)
        self.child_process.join(3)
        # Ensure the process is terminated
        self.assertIsNotNone(self.child_process)

    def remote_do(self, action):
        self.qsend.put(action)
        out = self.qrecv.get()
        self.assertNotIsInstance(out, BaseException)
        return out

    @staticmethod
    def do_terminate():
        return False, None


def remove_env(name):
    try:
        del os.environ[name]
    except KeyError:
        return False
    else:
        return True


def check_lib_lookup(qout, qin):
    status = True
    while status:
        try:
            action = qin.get()
        except Exception as e:
            qout.put(e)
            status = False
        else:
            try:
                with warnings.catch_warnings(record=True) as w:
                    warnings.simplefilter("always", NumbaWarning)
                    status, result = action()
                qout.put(result + (w,))
            except Exception as e:
                qout.put(e)
                status = False


@skip_on_cudasim('Library detection unsupported in the simulator')
@unittest.skipUnless(has_mp_get_context, 'mp.get_context not available')
@skip_unless_conda_cudatoolkit('test assumes conda installed cudatoolkit')
class TestLibDeviceLookUp(LibraryLookupBase):
    def test_libdevice_path_decision(self):
        # Check that the default is using conda environment
        by, info, warns = self.remote_do(self.do_clear_envs)
        if has_cuda:
            self.assertEqual(by, 'Conda environment')
        else:
            self.assertEqual(by, "<unknown>")
            self.assertIsNone(info)
        self.assertFalse(warns)
        # Check that CUDA_HOME works by removing conda-env
        by, info, warns = self.remote_do(self.do_set_cuda_home)
        self.assertEqual(by, 'CUDA_HOME')
        self.assertEqual(info, os.path.join('mycudahome', 'nvvm', 'libdevice'))
        self.assertFalse(warns)

        if get_system_ctk() is None:
            # Fake remove conda environment so no cudatoolkit is available
            by, info, warns = self.remote_do(self.do_clear_envs)
            self.assertEqual(by, '<unknown>')
            self.assertIsNone(info)
            self.assertFalse(warns)
        else:
            # Use system available cudatoolkit
            by, info, warns = self.remote_do(self.do_clear_envs)
            self.assertEqual(by, 'System')
            self.assertFalse(warns)

    @staticmethod
    def do_clear_envs():
        remove_env('CUDA_HOME')
        remove_env('CUDA_PATH')
        return True, _get_libdevice_path_decision()

    @staticmethod
    def do_set_cuda_home():
        os.environ['CUDA_HOME'] = os.path.join('mycudahome')
        _fake_non_conda_env()
        return True, _get_libdevice_path_decision()


@skip_on_cudasim('Library detection unsupported in the simulator')
@unittest.skipUnless(has_mp_get_context, 'mp.get_context not available')
@skip_unless_conda_cudatoolkit('test assumes conda installed cudatoolkit')
class TestNvvmLookUp(LibraryLookupBase):
    def test_nvvm_path_decision(self):
        # Check that the default is using conda environment
        by, info, warns = self.remote_do(self.do_clear_envs)
        if has_cuda:
            self.assertEqual(by, 'Conda environment')
        else:
            self.assertEqual(by, "<unknown>")
            self.assertIsNone(info)
        self.assertFalse(warns)
        # Check that CUDA_HOME works by removing conda-env
        by, info, warns = self.remote_do(self.do_set_cuda_home)
        self.assertEqual(by, 'CUDA_HOME')
        self.assertFalse(warns)
        if IS_WIN32:
            self.assertEqual(info, os.path.join('mycudahome', 'nvvm', 'bin'))
        elif IS_OSX:
            self.assertEqual(info, os.path.join('mycudahome', 'nvvm', 'lib'))
        else:
            self.assertEqual(info, os.path.join('mycudahome', 'nvvm', 'lib64'))

        if get_system_ctk() is None:
            # Fake remove conda environment so no cudatoolkit is available
            by, info, warns = self.remote_do(self.do_clear_envs)
            self.assertEqual(by, '<unknown>')
            self.assertIsNone(info)
            self.assertFalse(warns)
        else:
            # Use system available cudatoolkit
            by, info, warns = self.remote_do(self.do_clear_envs)
            self.assertEqual(by, 'System')
            self.assertFalse(warns)

    @staticmethod
    def do_clear_envs():
        remove_env('CUDA_HOME')
        remove_env('CUDA_PATH')
        return True, _get_nvvm_path_decision()

    @staticmethod
    def do_set_cuda_home():
        os.environ['CUDA_HOME'] = os.path.join('mycudahome')
        _fake_non_conda_env()
        return True, _get_nvvm_path_decision()


@skip_on_cudasim('Library detection unsupported in the simulator')
@unittest.skipUnless(has_mp_get_context, 'mp.get_context not available')
@skip_unless_conda_cudatoolkit('test assumes conda installed cudatoolkit')
class TestCudaLibLookUp(LibraryLookupBase):
    def test_cudalib_path_decision(self):
        # Check that the default is using conda environment
        by, info, warns = self.remote_do(self.do_clear_envs)
        if has_cuda:
            self.assertEqual(by, 'Conda environment')
        else:
            self.assertEqual(by, "<unknown>")
            self.assertIsNone(info)
        self.assertFalse(warns)

        # Check that CUDA_HOME works by removing conda-env
        self.remote_do(self.do_clear_envs)
        by, info, warns = self.remote_do(self.do_set_cuda_home)
        self.assertEqual(by, 'CUDA_HOME')
        self.assertFalse(warns)
        if IS_WIN32:
            self.assertEqual(info, os.path.join('mycudahome', 'bin'))
        elif IS_OSX:
            self.assertEqual(info, os.path.join('mycudahome', 'lib'))
        else:
            self.assertEqual(info, os.path.join('mycudahome', 'lib64'))
        if get_system_ctk() is None:
            # Fake remove conda environment so no cudatoolkit is available
            by, info, warns = self.remote_do(self.do_clear_envs)
            self.assertEqual(by, "<unknown>")
            self.assertIsNone(info)
            self.assertFalse(warns)
        else:
            # Use system available cudatoolkit
            by, info, warns = self.remote_do(self.do_clear_envs)
            self.assertEqual(by, 'System')
            self.assertFalse(warns)

    @staticmethod
    def do_clear_envs():
        remove_env('CUDA_HOME')
        remove_env('CUDA_PATH')
        return True, _get_cudalib_dir_path_decision()

    @staticmethod
    def do_set_cuda_home():
        os.environ['CUDA_HOME'] = os.path.join('mycudahome')
        _fake_non_conda_env()
        return True, _get_cudalib_dir_path_decision()


def _fake_non_conda_env():
    """
    Monkeypatch sys.prefix to hide the fact we are in a conda-env
    """
    sys.prefix = ''


if __name__ == '__main__':
    unittest.main()
