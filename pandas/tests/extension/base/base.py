import pandas.util.testing as tm


class BaseExtensionTests:
    @staticmethod
    def assert_series_equal(left, right, **kwargs):
        tm.assert_series_equal(left, right, **kwargs)

    @staticmethod
    def assert_frame_equal(left, right, **kwargs):
        tm.assert_frame_equal(left, right, **kwargs)
