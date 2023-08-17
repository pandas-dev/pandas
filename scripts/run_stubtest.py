import os
from pathlib import Path
import sys
import tempfile
import warnings

from mypy import stubtest

import pandas as pd

pd_version = getattr(pd, "__version__", "")

# fail early if pandas is not installed
if not pd_version:
    # fail on the CI, soft fail during local development
    warnings.warn("You need to install the development version of pandas")
    if pd.compat.is_ci_environment():
        sys.exit(1)
    else:
        sys.exit(0)

# GH 48260
if "dev" not in pd_version:
    warnings.warn(
        f"stubtest may fail as {pd_version} is not a dev version. "
        f"Please install a pandas dev version or see https://pandas.pydata.org/"
        f"pandas-docs/stable/development/contributing_codebase.html"
        f"#validating-type-hints on how to skip the stubtest"
    )


_ALLOWLIST = [  # should be empty
    # TODO (child classes implement these methods)
    "pandas._libs.hashtable.HashTable.__contains__",
    "pandas._libs.hashtable.HashTable.__len__",
    "pandas._libs.hashtable.HashTable.factorize",
    "pandas._libs.hashtable.HashTable.get_item",
    "pandas._libs.hashtable.HashTable.get_labels",
    "pandas._libs.hashtable.HashTable.get_na",
    "pandas._libs.hashtable.HashTable.get_state",
    "pandas._libs.hashtable.HashTable.lookup",
    "pandas._libs.hashtable.HashTable.map_locations",
    "pandas._libs.hashtable.HashTable.set_item",
    "pandas._libs.hashtable.HashTable.set_na",
    "pandas._libs.hashtable.HashTable.sizeof",
    "pandas._libs.hashtable.HashTable.unique",
    # stubtest might be too sensitive
    "pandas._libs.lib.NoDefault",
    "pandas._libs.lib._NoDefault.no_default",
    # internal type alias (should probably be private)
    "pandas._libs.lib.ndarray_obj_2d",
    # runtime argument "owner" has a default value but stub argument does not
    "pandas._libs.properties.AxisProperty.__get__",
    "pandas._libs.properties.cache_readonly.deleter",
    "pandas._libs.properties.cache_readonly.getter",
    "pandas._libs.properties.cache_readonly.setter",
    # TODO (child classes implement these methods)
    "pandas._libs.sparse.SparseIndex.__init__",
    "pandas._libs.sparse.SparseIndex.equals",
    "pandas._libs.sparse.SparseIndex.indices",
    "pandas._libs.sparse.SparseIndex.intersect",
    "pandas._libs.sparse.SparseIndex.lookup",
    "pandas._libs.sparse.SparseIndex.lookup_array",
    "pandas._libs.sparse.SparseIndex.make_union",
    "pandas._libs.sparse.SparseIndex.nbytes",
    "pandas._libs.sparse.SparseIndex.ngaps",
    "pandas._libs.sparse.SparseIndex.to_block_index",
    "pandas._libs.sparse.SparseIndex.to_int_index",
    # TODO (decorator changes argument names)
    "pandas._libs.tslibs.offsets.BaseOffset._apply_array",
    "pandas._libs.tslibs.offsets.BusinessHour.rollback",
    "pandas._libs.tslibs.offsets.BusinessHour.rollforward ",
    # type alias
    "pandas._libs.tslibs.timedeltas.UnitChoices",
]

if __name__ == "__main__":
    # find pyi files
    root = Path.cwd()
    pyi_modules = [
        str(pyi.relative_to(root).with_suffix("")).replace(os.sep, ".")
        for pyi in root.glob("pandas/**/*.pyi")
    ]

    # create allowlist
    with tempfile.NamedTemporaryFile(mode="w+t") as allow:
        allow.write("\n".join(_ALLOWLIST))
        allow.flush()

        args = pyi_modules + [
            "--ignore-missing-stub",
            "--concise",
            "--mypy-config-file",
            "pyproject.toml",
            "--allowlist",
            allow.name,
        ]
        sys.exit(stubtest.test_stubs(stubtest.parse_options(args)))
