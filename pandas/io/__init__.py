from typing import TYPE_CHECKING

if TYPE_CHECKING:
    # import modules that have public classes/functions
    from pandas.io import formats
    from pandas.io import json
    from pandas.io import stata

    # and mark only those modules as public
    __all__ = ["formats", "json", "stata"]
