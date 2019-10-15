from typing import Any

import pandas.util.testing as tm


class BaseExtensionTests:

    assert_equal = staticmethod(tm.assert_equal)
    assert_series_equal = staticmethod(tm.assert_series_equal)  # type: Any
    assert_frame_equal = staticmethod(tm.assert_frame_equal)  # type: Any
    assert_extension_array_equal = staticmethod(tm.assert_extension_array_equal)
