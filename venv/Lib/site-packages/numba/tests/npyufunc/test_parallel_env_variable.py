from numba.np.ufunc.parallel import get_thread_count
from os import environ as env
from numba.core import config
import unittest


class TestParallelEnvVariable(unittest.TestCase):
    """
    Tests environment variables related to the underlying "parallel"
    functions for npyufuncs.
    """

    _numba_parallel_test_ = False

    def test_num_threads_variable(self):
        """
        Tests the NUMBA_NUM_THREADS env variable behaves as expected.
        """
        key = 'NUMBA_NUM_THREADS'
        current = str(getattr(env, key, config.NUMBA_NUM_THREADS))
        threads = "3154"
        env[key] = threads
        try:
            config.reload_config()
        except RuntimeError as e:
            # This test should fail if threads have already been launched
            self.assertIn("Cannot set NUMBA_NUM_THREADS", e.args[0])
        else:
            self.assertEqual(threads, str(get_thread_count()))
            self.assertEqual(threads, str(config.NUMBA_NUM_THREADS))
        finally:
            # reset the env variable/set to default. Should not fail even if
            # threads are launched because the value is the same.
            env[key] = current
            config.reload_config()

if __name__ == '__main__':
    unittest.main()
