import os
import tempfile
from textwrap import dedent
import unittest
from unittest import mock
from numba.tests.support import (TestCase, temp_directory, override_env_config,
                                 run_in_subprocess)
from numba.core import config

try:
    import yaml
    _HAVE_YAML = True
except ImportError:
    _HAVE_YAML = False

_skip_msg = "pyyaml needed for configuration file tests"
needs_yaml = unittest.skipIf(not _HAVE_YAML, _skip_msg)


@needs_yaml
class TestConfig(TestCase):

    # Disable parallel testing due to envvars modification
    _numba_parallel_test_ = False

    def setUp(self):
        # use support.temp_directory, it can do the clean up
        self.tmppath = temp_directory('config_tmp')
        self.maxDiff = 2500
        super(TestConfig, self).setUp()

    def mock_cfg_location(self):
        """
        Creates a mock launch location.
        Returns the location path.
        """
        return tempfile.mkdtemp(dir=self.tmppath)

    def inject_mock_cfg(self, location, cfg):
        """
        Injects a mock configuration at 'location'
        """
        tmpcfg = os.path.join(location, config._config_fname)
        with open(tmpcfg, 'wt') as f:
            yaml.dump(cfg, f, default_flow_style=False)

    def get_settings(self):
        """
        Gets the current numba config settings
        """
        store = dict()
        for x in dir(config):
            if x.isupper():
                store[x] = getattr(config, x)
        return store

    def create_config_effect(self, cfg):
        """
        Returns a config "original" from a location with no config file
        and then the impact of applying the supplied cfg dictionary as
        a config file at a location in the returned "current".
        """

        # store original cwd
        original_cwd = os.getcwd()

        # create mock launch location
        launch_dir = self.mock_cfg_location()

        # switch cwd to the mock launch location, get and store settings
        os.chdir(launch_dir)
        # use override to ensure that the config is zero'd out with respect
        # to any existing settings
        with override_env_config('_', '_'):
            original = self.get_settings()

        # inject new config into a file in the mock launch location
        self.inject_mock_cfg(launch_dir, cfg)

        try:
            # override something but don't change the value, this is to refresh
            # the config and make sure the injected config file is read
            with override_env_config('_', '_'):
                current = self.get_settings()
        finally:
            # switch back to original dir with no new config
            os.chdir(original_cwd)
        return original, current

    def test_config(self):
        # ensure a non empty settings file does impact config and that the
        # case of the key makes no difference
        key = 'COLOR_SCHEME'
        for case in [str.upper, str.lower]:
            orig, curr = self.create_config_effect({case(key): 'light_bg'})
            self.assertTrue(orig != curr)
            self.assertTrue(orig[key] != curr[key])
            self.assertEqual(curr[key], 'light_bg')
            # check that just the color scheme is the cause of difference
            orig.pop(key)
            curr.pop(key)
            self.assertEqual(orig, curr)

    def test_empty_config(self):
        # ensure an empty settings file does not impact config
        orig, curr = self.create_config_effect({})
        self.assertEqual(orig, curr)

    def test_illegal_error_style_handling(self):
        # ensure that illegal error styles are ignored
        new_env = os.environ.copy()
        new_env['NUMBA_CAPTURED_ERRORS'] = 'not_a_known_style'
        source_compiled = "the source compiled"
        code = ("from numba import njit\n@njit\ndef foo():\n\t"
                f"print('{source_compiled}')\nfoo()")
        out, err = run_in_subprocess(dedent(code), env=new_env)
        expected = ("Environment variable \'NUMBA_CAPTURED_ERRORS\' is defined "
                    "but its associated value \'not_a_known_style\' could not "
                    "be parsed.")
        err_msg = err.decode('utf-8')
        self.assertIn(expected, err_msg)
        ex_expected = ("Invalid style in NUMBA_CAPTURED_ERRORS: "
                       "not_a_known_style")
        self.assertIn(ex_expected, err_msg)
        self.assertIn(source_compiled, out.decode('utf-8'))

    def test_default_error_style_handling(self):
        # ensure that the default is new_style
        new_env = os.environ.copy()
        new_env['NUMBA_CAPTURED_ERRORS'] = 'default'
        code = ("from numba.core import config\n"
                "print('---->', config.CAPTURED_ERRORS)\n"
                "assert config.CAPTURED_ERRORS == 'new_style'")
        out, err = run_in_subprocess(dedent(code), env=new_env)
        err_msg = err.decode('utf-8')
        out_msg = out.decode('utf-8')
        ex_expected = "----> new_style"
        self.assertIn(ex_expected, out_msg, msg=err_msg)

    @unittest.skipUnless(config.ENABLE_AVX,
                         "test expects NUMBA_ENABLE_AVX==True")
    def test_nocona_disables_avx(self):
        # test with nocona
        new_env = os.environ.copy()
        new_env.pop('NUMBA_ENABLE_AVX', None)  # clear NUMBA_ENABLE_AVX

        new_env['NUMBA_CPU_NAME'] = 'nocona'
        code = ("from numba.core import config\n"
                "print('---->', bool(config.ENABLE_AVX))\n"
                "assert not config.ENABLE_AVX")
        out, err = run_in_subprocess(dedent(code), env=new_env)
        err_msg = err.decode('utf-8')
        out_msg = out.decode('utf-8')
        ex_expected = "----> False"
        self.assertIn(ex_expected, out_msg, msg=err_msg)

        # test with skylake-avx512
        new_env['NUMBA_CPU_NAME'] = 'skylake-avx512'
        code = ("from numba.core import config\n"
                "print('---->', bool(config.ENABLE_AVX))\n"
                "assert config.ENABLE_AVX")
        out, err = run_in_subprocess(dedent(code), env=new_env)
        err_msg = err.decode('utf-8')
        out_msg = out.decode('utf-8')
        ex_expected = "----> True"
        self.assertIn(ex_expected, out_msg, msg=err_msg)


class TestNumbaOptLevel(TestCase):
    # Tests that the setting of NUMBA_OPT influences the "cheap" module pass.
    # Spot checks NUMBA_OPT={'max', '3', '0'}

    def check(self, expected, opt_value, raw_value):
        # local imports for state-safety
        from numba import config, njit

        # check opt value and its raw_value
        self.assertEqual(config.OPT, opt_value)
        self.assertEqual(config.OPT._raw_value, raw_value)

        # Patch the CPUCodegen to make capture calls to the
        # `_module_pass_manager` through a `side_effect` function that asserts
        # that the kwargs being passed are as expected per the "NUMBA_OPT"
        # level. The `side_effect` function immediately raises with a knwon
        # message to abort further stages compilation once the check is
        # complete.
        from numba.core.codegen import CPUCodegen
        side_effect_message = "expected side effect"

        def side_effect(*args, **kwargs):
            self.assertEqual(kwargs, expected)
            raise RuntimeError(side_effect_message)

        with mock.patch.object(CPUCodegen, '_module_pass_manager',
                               side_effect=side_effect):
            with self.assertRaises(RuntimeError) as raises:
                njit(lambda : ...)()

            self.assertIn(side_effect_message, str(raises.exception))

    @TestCase.run_test_in_subprocess(envvars={'NUMBA_OPT': 'max'})
    def test_opt_max(self):
        # NUMBA_OPT='max' should set opt to 3 and enable loop_vectorize
        expected = {'loop_vectorize': True,
                    'slp_vectorize': False,
                    'opt': 3,
                    'cost': 'cheap'}
        self.check(expected, 3, 'max')

    @TestCase.run_test_in_subprocess(envvars={'NUMBA_OPT': '3'})
    def test_opt_3(self):
        # NUMBA_OPT='3' should not impact opt or loop_vectorize
        expected = {'loop_vectorize': False,
                    'slp_vectorize': False,
                    'opt': 0,
                    'cost': 'cheap'}
        self.check(expected, 3, 3)

    @TestCase.run_test_in_subprocess(envvars={'NUMBA_OPT': '0'})
    def test_opt_0(self):
        # NUMBA_OPT='0' should not impact opt or loop_vectorize
        expected = {'loop_vectorize': False,
                    'slp_vectorize': False,
                    'opt': 0,
                    'cost': 'cheap'}
        self.check(expected, 0, 0)

    @TestCase.run_test_in_subprocess()
    def test_opt_default(self):
        # NUMBA_OPT is not set, the default should not impact opt or
        # loop_vectorize
        expected = {'loop_vectorize': False,
                    'slp_vectorize': False,
                    'opt': 0,
                    'cost': 'cheap'}
        self.check(expected, 3, 3)

    @TestCase.run_test_in_subprocess(envvars={'NUMBA_OPT': 'invalid'})
    def test_opt_invalid(self):
        # NUMBA_OPT='invalid' should just proceed as default case
        expected = {'loop_vectorize': False,
                    'slp_vectorize': False,
                    'opt': 0,
                    'cost': 'cheap'}
        self.check(expected, 3, 3)


if __name__ == '__main__':
    unittest.main()
