import os
from pathlib import Path
import sys
import tempfile
import warnings

from mypy import stubtest

import pandas as pd

# fail early if pandas is not installed
if "dev" not in getattr(pd, "__version__", ""):
    # 'succeed' but print warning
    warnings.warn("You need to install the development version of pandas")
    sys.exit(0)


_ALLOWLIST = [  # should be empty
    # TODO
    "pandas._libs.hashtable.Complex128Vector.__init__",
    "pandas._libs.hashtable.Complex64Vector.__init__",
    "pandas._libs.hashtable.Float32Vector.__init__",
    "pandas._libs.hashtable.Float64Vector.__init__",
    "pandas._libs.hashtable.HashTable.__contains__",
    "pandas._libs.hashtable.HashTable.__len__",
    "pandas._libs.hashtable.HashTable.factorize",
    "pandas._libs.hashtable.HashTable.get_item",
    "pandas._libs.hashtable.HashTable.get_labels",
    "pandas._libs.hashtable.HashTable.get_state",
    "pandas._libs.hashtable.HashTable.lookup",
    "pandas._libs.hashtable.HashTable.map_locations",
    "pandas._libs.hashtable.HashTable.set_item",
    "pandas._libs.hashtable.HashTable.sizeof",
    "pandas._libs.hashtable.HashTable.unique",
    "pandas._libs.hashtable.Int16Vector.__init__",
    "pandas._libs.hashtable.Int32Vector.__init__",
    "pandas._libs.hashtable.Int64Vector.__init__",
    "pandas._libs.hashtable.Int8Vector.__init__",
    "pandas._libs.hashtable.ObjectVector.__init__",
    "pandas._libs.hashtable.StringVector.__init__",
    "pandas._libs.hashtable.UInt16Vector.__init__",
    "pandas._libs.hashtable.UInt32Vector.__init__",
    "pandas._libs.hashtable.UInt64Vector.__init__",
    "pandas._libs.hashtable.UInt8Vector.__init__",
    # stubtest might be too sensitive
    "pandas._libs.lib.NoDefault",
    "pandas._libs.lib._NoDefault.no_default",
    # internal type alias (should probably be private)
    "pandas._libs.lib.ndarray_obj_2d",
    # workaround for mypy (cache_readonly = property)
    "pandas._libs.properties.cache_readonly.__get__",
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
