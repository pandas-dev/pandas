###################################################################
#  Numexpr - Fast numerical array expression evaluator for NumPy.
#
#      License: MIT
#      Author:  See AUTHORS.txt
#
#  See LICENSE.txt and LICENSES/*.txt for details about copyright and
#  rights to use.
####################################################################

import logging
log = logging.getLogger(__name__)

import os
import subprocess

from numexpr.interpreter import _set_num_threads, _get_num_threads, MAX_THREADS
from numexpr import use_vml
from . import version

if use_vml:
    from numexpr.interpreter import (
        _get_vml_version, _set_vml_accuracy_mode, _set_vml_num_threads,
        _get_vml_num_threads)


def get_vml_version():
    """
    Get the VML/MKL library version.
    """
    if use_vml:
        return _get_vml_version()
    else:
        return None


def set_vml_accuracy_mode(mode):
    """
    Set the accuracy mode for VML operations.

    The `mode` parameter can take the values:
    - 'high': high accuracy mode (HA), <1 least significant bit
    - 'low': low accuracy mode (LA), typically 1-2 least significant bits
    - 'fast': enhanced performance mode (EP)
    - None: mode settings are ignored

    This call is equivalent to the `vmlSetMode()` in the VML library.
    See:

    http://www.intel.com/software/products/mkl/docs/webhelp/vml/vml_DataTypesAccuracyModes.html

    for more info on the accuracy modes.

    Returns old accuracy settings.
    """
    if use_vml:
        acc_dict = {None: 0, 'low': 1, 'high': 2, 'fast': 3}
        acc_reverse_dict = {1: 'low', 2: 'high', 3: 'fast'}
        if mode not in list(acc_dict.keys()):
            raise ValueError(
                "mode argument must be one of: None, 'high', 'low', 'fast'")
        retval = _set_vml_accuracy_mode(acc_dict.get(mode, 0))
        return acc_reverse_dict.get(retval)
    else:
        return None


def set_vml_num_threads(nthreads):
    """
    Suggests a maximum number of threads to be used in VML operations.

    This function is equivalent to the call
    `mkl_domain_set_num_threads(nthreads, MKL_DOMAIN_VML)` in the MKL
    library.  See:

    http://www.intel.com/software/products/mkl/docs/webhelp/support/functn_mkl_domain_set_num_threads.html

    for more info about it.
    """
    if use_vml:
        _set_vml_num_threads(nthreads)
    pass

def get_vml_num_threads():
    """
    Gets the maximum number of threads to be used in VML operations.

    This function is equivalent to the call
    `mkl_domain_get_max_threads (MKL_DOMAIN_VML)` in the MKL
    library.  See:

    http://software.intel.com/en-us/node/522118

    for more info about it.
    """
    if use_vml:
        return _get_vml_num_threads()
    return None

def set_num_threads(nthreads):
    """
    Sets a number of threads to be used in operations.

    DEPRECATED: returns the previous setting for the number of threads.

    During initialization time NumExpr sets this number to the number
    of detected cores in the system (see `detect_number_of_cores()`).
    """
    old_nthreads = _set_num_threads(nthreads)
    return old_nthreads

def get_num_threads():
    """
    Gets the number of threads currently in use for operations.
    """
    return _get_num_threads()

def _init_num_threads():
    """
    Detects the environment variable 'NUMEXPR_MAX_THREADS' to set the threadpool 
    size, and if necessary the slightly redundant 'NUMEXPR_NUM_THREADS' or 
    'OMP_NUM_THREADS' env vars to set the initial number of threads used by 
    the virtual machine.
    """
    # Any platform-specific short-circuits
    if 'sparc' in version.platform_machine:
        log.warning('The number of threads have been set to 1 because problems related '
                  'to threading have been reported on some sparc machine. '
                  'The number of threads can be changed using the "set_num_threads" '
                  'function.')
        set_num_threads(1)
        return 1

    env_configured = False
    n_cores = detect_number_of_cores()
    if ('NUMEXPR_MAX_THREADS' in os.environ and os.environ['NUMEXPR_MAX_THREADS'] != '' or
        'OMP_NUM_THREADS' in os.environ and os.environ['OMP_NUM_THREADS'] != ''):
        # The user has configured NumExpr in the expected way, so suppress logs.
        env_configured = True
        n_cores = MAX_THREADS
    else:
        # The use has not set 'NUMEXPR_MAX_THREADS', so likely they have not 
        # configured NumExpr as desired, so we emit info logs.
        if n_cores > MAX_THREADS:
            log.info('Note: detected %d virtual cores but NumExpr set to maximum of %d, check "NUMEXPR_MAX_THREADS" environment variable.'%(n_cores, MAX_THREADS))
        if n_cores > 16:
            # Back in 2019, 8 threads would be considered safe for performance. We are in 2024 now, so adjusting.
            log.info('Note: NumExpr detected %d cores but "NUMEXPR_MAX_THREADS" not set, so enforcing safe limit of 16.'%n_cores)
            n_cores = 16

    # Now we check for 'NUMEXPR_NUM_THREADS' or 'OMP_NUM_THREADS' to set the 
    # actual number of threads used.
    if 'NUMEXPR_NUM_THREADS' in os.environ and os.environ['NUMEXPR_NUM_THREADS'] != '':
        requested_threads = int(os.environ['NUMEXPR_NUM_THREADS'])
    elif 'OMP_NUM_THREADS' in os.environ and os.environ['OMP_NUM_THREADS'] != '':
        # Empty string is commonly used to unset the variable
        requested_threads = int(os.environ['OMP_NUM_THREADS'])
    else:
        requested_threads = n_cores
        if not env_configured:
            log.info('NumExpr defaulting to %d threads.'%n_cores)

    # The C-extension function performs its own checks against `MAX_THREADS`
    set_num_threads(requested_threads)
    return requested_threads

    
def detect_number_of_cores():
    """
    Detects the number of cores on a system. Cribbed from pp.
    """
    # Linux, Unix and MacOS:
    if hasattr(os, "sysconf"):
        if "SC_NPROCESSORS_ONLN" in os.sysconf_names:
            # Linux & Unix:
            ncpus = os.sysconf("SC_NPROCESSORS_ONLN")
            if isinstance(ncpus, int) and ncpus > 0:
                return ncpus
        else:  # OSX:
            return int(subprocess.check_output(["sysctl", "-n", "hw.ncpu"]))
    # Windows:
    try:
        ncpus = int(os.environ.get("NUMBER_OF_PROCESSORS", ""))
        if ncpus > 0:
            return ncpus
    except ValueError:
        pass
    return 1  # Default


def detect_number_of_threads():
    """
    DEPRECATED: use `_init_num_threads` instead.
    If this is modified, please update the note in: https://github.com/pydata/numexpr/wiki/Numexpr-Users-Guide
    """
    log.warning('Deprecated, use `init_num_threads` instead.')
    try:
        nthreads = int(os.environ.get('NUMEXPR_NUM_THREADS', ''))
    except ValueError:
        try:
            nthreads = int(os.environ.get('OMP_NUM_THREADS', ''))
        except ValueError:
            nthreads = detect_number_of_cores()

    # Check that we don't surpass the MAX_THREADS in interpreter.cpp
    if nthreads > MAX_THREADS:
        nthreads = MAX_THREADS
    return nthreads


class CacheDict(dict):
    """
    A dictionary that prevents itself from growing too much.
    """

    def __init__(self, maxentries):
        self.maxentries = maxentries
        super(CacheDict, self).__init__(self)

    def __setitem__(self, key, value):
        # Protection against growing the cache too much
        if len(self) > self.maxentries:
            # Remove a 10% of (arbitrary) elements from the cache
            entries_to_remove = self.maxentries // 10
            for k in list(self.keys())[:entries_to_remove]:
                super(CacheDict, self).__delitem__(k)
        super(CacheDict, self).__setitem__(key, value)

