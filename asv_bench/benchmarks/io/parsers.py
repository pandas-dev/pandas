try:
    from pandas._libs.tslibs.parsing import _does_string_look_like_datetime
except ImportError:
    # Avoid whole benchmark suite import failure on asv (currently 0.4)
    pass


class DoesStringLookLikeDatetime:
    params = (["2Q2005", "0.0", "10000"],)
    param_names = ["value"]

    def setup(self, value):
        self.objects = [value] * 1000000

    def time_check_datetimes(self, value):
        for obj in self.objects:
            _does_string_look_like_datetime(obj)
