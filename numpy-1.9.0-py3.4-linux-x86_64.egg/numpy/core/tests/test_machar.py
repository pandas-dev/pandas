from __future__ import division, absolute_import, print_function

from numpy.testing import *

from numpy.core.machar import MachAr
import numpy.core.numerictypes as ntypes
from numpy import errstate, array

class TestMachAr(TestCase):
    def _run_machar_highprec(self):
        # Instanciate MachAr instance with high enough precision to cause
        # underflow
        try:
            hiprec = ntypes.float96
            machar = MachAr(lambda v:array([v], hiprec))
        except AttributeError:
            "Skipping test: no nyptes.float96 available on this platform."

    def test_underlow(self):
        """Regression testing for #759: instanciating MachAr for dtype =
        np.float96 raises spurious warning."""
        with errstate(all='raise'):
            try:
                self._run_machar_highprec()
            except FloatingPointError as e:
                self.fail("Caught %s exception, should not have been raised." % e)


if __name__ == "__main__":
    run_module_suite()
