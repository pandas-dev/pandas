# Licensed under a 3-clause BSD style license - see LICENSE.rst

class TimerawSuite:
    def timeraw_fresh(self):
        # The test interpreter should not be polluted with anything.
        return """
        import sys
        assert 'asv' not in sys.modules
        assert 'inspect' not in sys.modules
        """

    def timeraw_setup(self):
        # Setup code
        return "a+1", "a=3"

    def timeraw_timeout(self):
        return """
        # Inside the grandchild.
        while True:
            pass
        """
    timeraw_timeout.timeout = 0.1

    def timeraw_count(self):
        # Using number other than one or a fixed repeat should work.
        return """
        import sys
        sys.stderr.write('0')
        """
    timeraw_count.repeat = 7
    timeraw_count.number = 3
    timeraw_count.warmup_time = 0
