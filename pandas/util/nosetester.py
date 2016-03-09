"""
Nose test running.

This module implements ``test()`` function for pandas modules.

"""
from __future__ import division, absolute_import, print_function

import os
import sys
import warnings
from pandas.compat import string_types
from numpy.testing import nosetester


def get_package_name(filepath):
    """
    Given a path where a package is installed, determine its name.

    Parameters
    ----------
    filepath : str
        Path to a file. If the determination fails, "pandas" is returned.

    Examples
    --------
    >>> pandas.util.nosetester.get_package_name('nonsense')
    'pandas'

    """

    pkg_name = []
    while 'site-packages' in filepath or 'dist-packages' in filepath:
        filepath, p2 = os.path.split(filepath)
        if p2 in ('site-packages', 'dist-packages'):
            break
        pkg_name.append(p2)

    # if package name determination failed, just default to pandas
    if not pkg_name:
        return "pandas"

    # otherwise, reverse to get correct order and return
    pkg_name.reverse()

    # don't include the outer egg directory
    if pkg_name[0].endswith('.egg'):
        pkg_name.pop(0)

    return '.'.join(pkg_name)

import_nose = nosetester.import_nose
run_module_suite = nosetester.run_module_suite


class NoseTester(nosetester.NoseTester):
    """
    Nose test runner.

    This class is made available as pandas.util.nosetester.NoseTester, and
    a test function is typically added to a package's __init__.py like so::

      from numpy.testing import Tester
      test = Tester().test

    Calling this test function finds and runs all tests associated with the
    package and all its sub-packages.

    Attributes
    ----------
    package_path : str
        Full path to the package to test.
    package_name : str
        Name of the package to test.

    Parameters
    ----------
    package : module, str or None, optional
        The package to test. If a string, this should be the full path to
        the package. If None (default), `package` is set to the module from
        which `NoseTester` is initialized.
    raise_warnings : None, str or sequence of warnings, optional
        This specifies which warnings to configure as 'raise' instead
        of 'warn' during the test execution.  Valid strings are:

          - "develop" : equals ``(DeprecationWarning, RuntimeWarning)``
          - "release" : equals ``()``, don't raise on any warnings.

        See Notes for more details.

    Notes
    -----
    The default for `raise_warnings` is
    ``(DeprecationWarning, RuntimeWarning)`` for development versions of
    pandas, and ``()`` for released versions.  The purpose of this switching
    behavior is to catch as many warnings as possible during development, but
    not give problems for packaging of released versions.

    """
    excludes = []

    def _show_system_info(self):
        nose = import_nose()

        import pandas
        print("pandas version %s" % pandas.__version__)
        import numpy
        print("numpy version %s" % numpy.__version__)
        pddir = os.path.dirname(pandas.__file__)
        print("pandas is installed in %s" % pddir)

        pyversion = sys.version.replace('\n', '')
        print("Python version %s" % pyversion)
        print("nose version %d.%d.%d" % nose.__versioninfo__)

    def _get_custom_doctester(self):
        """ Return instantiated plugin for doctests

        Allows subclassing of this class to override doctester

        A return value of None means use the nose builtin doctest plugin
        """
        return None

    def _test_argv(self, label, verbose, extra_argv):
        """
        Generate argv for nosetest command

        Parameters
        ----------
        label : {'fast', 'full', '', attribute identifier}, optional
            see ``test`` docstring
        verbose : int, optional
            Verbosity value for test outputs, in the range 1-10. Default is 1.
        extra_argv : list, optional
            List with any extra arguments to pass to nosetests.

        Returns
        -------
        argv : list
            command line arguments that will be passed to nose
        """

        argv = [__file__, self.package_path]
        if label and label != 'full':
            if not isinstance(label, string_types):
                raise TypeError('Selection label should be a string')
            if label == 'fast':
                label = 'not slow and not network and not disabled'
            argv += ['-A', label]
        argv += ['--verbosity', str(verbose)]

        # When installing with setuptools, and also in some other cases, the
        # test_*.py files end up marked +x executable. Nose, by default, does
        # not run files marked with +x as they might be scripts. However, in
        # our case nose only looks for test_*.py files under the package
        # directory, which should be safe.
        argv += ['--exe']

        if extra_argv:
            argv += extra_argv
        return argv

    def test(self, label='fast', verbose=1, extra_argv=None,
             doctests=False, coverage=False, raise_warnings=None):
        """
        Run tests for module using nose.

        Parameters
        ----------
        label : {'fast', 'full', '', attribute identifier}, optional
            Identifies the tests to run. This can be a string to pass to
            the nosetests executable with the '-A' option, or one of several
            special values.  Special values are:

            * 'fast' - the default - which corresponds to the ``nosetests -A``
              option of 'not slow'.
            * 'full' - fast (as above) and slow tests as in the
              'no -A' option to nosetests - this is the same as ''.
            * None or '' - run all tests.
            * attribute_identifier - string passed directly to nosetests
              as '-A'.

        verbose : int, optional
            Verbosity value for test outputs, in the range 1-10. Default is 1.
        extra_argv : list, optional
            List with any extra arguments to pass to nosetests.
        doctests : bool, optional
            If True, run doctests in module. Default is False.
        coverage : bool, optional
            If True, report coverage of NumPy code. Default is False.
            (This requires the `coverage module
            <http://nedbatchelder.com/code/modules/coverage.html>`_).
        raise_warnings : str or sequence of warnings, optional
            This specifies which warnings to configure as 'raise' instead
            of 'warn' during the test execution.  Valid strings are:

            - 'develop' : equals ``(DeprecationWarning, RuntimeWarning)``
            - 'release' : equals ``()``, don't raise on any warnings.

        Returns
        -------
        result : object
            Returns the result of running the tests as a
            ``nose.result.TextTestResult`` object.

        """

        # cap verbosity at 3 because nose becomes *very* verbose beyond that
        verbose = min(verbose, 3)

        if doctests:
            print("Running unit tests and doctests for %s" % self.package_name)
        else:
            print("Running unit tests for %s" % self.package_name)

        self._show_system_info()

        # reset doctest state on every run
        import doctest
        doctest.master = None

        if raise_warnings is None:

            # default based on if we are released
            from pandas import __version__
            from distutils.version import StrictVersion
            try:
                StrictVersion(__version__)
                raise_warnings = 'release'
            except ValueError:
                raise_warnings = 'develop'

        _warn_opts = dict(develop=(DeprecationWarning, RuntimeWarning),
                          release=())
        if isinstance(raise_warnings, string_types):
            raise_warnings = _warn_opts[raise_warnings]

        with warnings.catch_warnings():

            if len(raise_warnings):

                # Reset the warning filters to the default state,
                # so that running the tests is more repeatable.
                warnings.resetwarnings()
                # Set all warnings to 'warn', this is because the default
                # 'once' has the bad property of possibly shadowing later
                # warnings.
                warnings.filterwarnings('always')
                # Force the requested warnings to raise
                for warningtype in raise_warnings:
                    warnings.filterwarnings('error', category=warningtype)
                # Filter out annoying import messages.
                warnings.filterwarnings("ignore", category=FutureWarning)

            from numpy.testing.noseclasses import NumpyTestProgram
            argv, plugins = self.prepare_test_args(
                label, verbose, extra_argv, doctests, coverage)
            t = NumpyTestProgram(argv=argv, exit=False, plugins=plugins)

        return t.result
