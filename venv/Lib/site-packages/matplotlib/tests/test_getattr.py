from importlib import import_module
from pkgutil import walk_packages
import sys
import warnings

import pytest

import matplotlib
from matplotlib.testing import is_ci_environment, subprocess_run_helper

# Get the names of all matplotlib submodules,
# except for the unit tests and private modules.
module_names = []
backend_module_names = []
for m in walk_packages(path=matplotlib.__path__, prefix=f'{matplotlib.__name__}.'):
    if m.name.startswith(__package__):
        continue
    if any(x.startswith('_') for x in m.name.split('.')):
        continue
    if 'backends.backend_' in m.name:
        backend_module_names.append(m.name)
    else:
        module_names.append(m.name)


def _test_getattr(module_name, use_pytest=True):
    """
    Test that __getattr__ methods raise AttributeError for unknown keys.
    See #20822, #20855.
    """
    try:
        module = import_module(module_name)
    except (ImportError, RuntimeError, OSError) as e:
        # Skip modules that cannot be imported due to missing dependencies
        if use_pytest:
            pytest.skip(f'Cannot import {module_name} due to {e}')
        else:
            print(f'SKIP: Cannot import {module_name} due to {e}')
            return

    key = 'THIS_SYMBOL_SHOULD_NOT_EXIST'
    if hasattr(module, key):
        delattr(module, key)


@pytest.mark.parametrize('module_name', module_names)
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
@pytest.mark.filterwarnings('ignore::ImportWarning')
def test_getattr(module_name):
    _test_getattr(module_name)


def _test_module_getattr():
    warnings.filterwarnings('ignore', category=DeprecationWarning)
    warnings.filterwarnings('ignore', category=ImportWarning)
    module_name = sys.argv[1]
    _test_getattr(module_name, use_pytest=False)


@pytest.mark.parametrize('module_name', backend_module_names)
def test_backend_getattr(module_name):
    proc = subprocess_run_helper(_test_module_getattr, module_name,
                                 timeout=120 if is_ci_environment() else 20)
    if 'SKIP: ' in proc.stdout:
        pytest.skip(proc.stdout.removeprefix('SKIP: '))
    print(proc.stdout)
