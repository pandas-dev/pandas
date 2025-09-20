import unittest
import pickle
import sys
import tempfile
from pathlib import Path


class TestUnpickleDeletedModule(unittest.TestCase):
    def test_loading_pickle_with_no_module(self):
        """Create a module that uses Numba, import a function from it.
        Then delete the module and pickle the function. The function
        should load from the pickle without a problem.

        Note - This is a simplified version of how Numba might be used
        on a distributed system using e.g. dask distributed. With the
        pickle being sent to the worker but not the original module.
        """

        # Source code for temporary module we will make
        source = "\n".join(
            [
                "from numba import vectorize",
                "@vectorize(['float64(float64)'])",
                "def inc1(x):",
                "    return x + 1",
            ]
        )

        # Create a temporary directory and add it to path.
        modname = "tmp_module"
        with tempfile.TemporaryDirectory() as tmp_dir:
            sys.path.append(tmp_dir)

            # Create tmp_module.py in there with our source code above.
            filename = Path(f"{tmp_dir}/{modname}.py")
            f = open(filename, "a")
            f.write(source)
            f.close()

            # Import the temporary module before file is deleted
            from tmp_module import inc1

        # Remove from imported libraries
        del sys.modules[modname]

        # Pickle function and assert that it loads correctly
        pkl = pickle.dumps(inc1)
        f = pickle.loads(pkl)
        self.assertEqual(f(2), 3)
