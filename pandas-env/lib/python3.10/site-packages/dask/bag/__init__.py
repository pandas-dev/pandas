from __future__ import annotations

try:
    from dask.bag.avro import read_avro
    from dask.bag.core import Bag, Item
    from dask.bag.core import bag_map as map
    from dask.bag.core import bag_range as range
    from dask.bag.core import bag_zip as zip
    from dask.bag.core import (
        concat,
        from_delayed,
        from_sequence,
        from_url,
        map_partitions,
        to_textfiles,
    )
    from dask.bag.text import read_text
    from dask.bag.utils import assert_eq
    from dask.base import compute
except ImportError as e:
    msg = (
        "Dask bag requirements are not installed.\n\n"
        "Please either conda or pip install as follows:\n\n"
        "  conda install dask               # either conda install\n"
        '  python -m pip install "dask[bag]" --upgrade  # or python -m pip install'
    )
    raise ImportError(str(e) + "\n\n" + msg) from e
