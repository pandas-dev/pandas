import pandas._testing as tm
from pandas.core.arrays import (
    ExtensionArray,
    ExtensionOpsMixin,
    ExtensionScalarOpsMixin,
)


def test_extension_ops_mixin_deprecated():
    # GH#37080 deprecated in favor of OpsMixin
    with tm.assert_produces_warning(FutureWarning):

        class MySubclass(ExtensionOpsMixin, ExtensionArray):
            pass

    with tm.assert_produces_warning(FutureWarning):

        class MyOtherSubclass(ExtensionScalarOpsMixin, ExtensionArray):
            pass
