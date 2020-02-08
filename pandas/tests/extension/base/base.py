import pandas._testing as tm


class BaseExtensionTests:
    # can't use staticmethod() as this confuses mypy
    def assert_equal(self, left, right, **kwargs):
        return tm.assert_equal(left, right, **kwargs)

    def assert_series_equal(self, left, right, *args, **kwargs):
        return tm.assert_series_equal(left, right, *args, **kwargs)

    def assert_frame_equal(self, left, right, *args, **kwargs):
        return tm.assert_frame_equal(left, right, *args, **kwargs)

    def assert_extension_array_equal(self, left, right, *args, **kwargs):
        return tm.assert_extension_array_equal(left, right, *args, **kwargs)
