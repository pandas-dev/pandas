from __future__ import annotations

import functools

import numpy as np
import toolz

from dask.dataframe.dask_expr._expr import Expr
from dask.dataframe.partitionquantiles import (
    create_merge_tree,
    dtype_info,
    merge_and_compress_summaries,
    percentiles_summary,
    process_val_weights,
)
from dask.tokenize import tokenize
from dask.utils import random_state_data


class RepartitionQuantiles(Expr):
    _parameters = ["frame", "input_npartitions", "upsample", "random_state"]
    _defaults = {"upsample": 1.0, "random_state": None}

    @functools.cached_property
    def _meta(self):
        return self.frame._meta

    @property
    def npartitions(self):
        return 1

    def _divisions(self):
        return 0.0, 1.0

    def __dask_postcompute__(self):
        return toolz.first, ()

    def _layer(self):
        import pandas as pd

        qs = np.linspace(0, 1, self.input_npartitions + 1)
        if self.random_state is None:
            random_state = int(tokenize(self.operands), 16) % np.iinfo(np.int32).max
        else:
            random_state = self.random_state
        state_data = random_state_data(self.frame.npartitions, random_state)

        keys = self.frame.__dask_keys__()
        dtype_dsk = {(self._name, 0, 0): (dtype_info, keys[0])}

        percentiles_dsk = {
            (self._name, 1, i): (
                percentiles_summary,
                key,
                self.frame.npartitions,
                self.input_npartitions,
                self.upsample,
                state,
            )
            for i, (state, key) in enumerate(zip(state_data, keys))
        }

        merge_dsk = create_merge_tree(
            merge_and_compress_summaries, sorted(percentiles_dsk), self._name, 2
        )
        if not merge_dsk:
            # Compress the data even if we only have one partition
            merge_dsk = {
                (self._name, 2, 0): (
                    merge_and_compress_summaries,
                    [list(percentiles_dsk)[0]],
                )
            }

        merged_key = max(merge_dsk)
        last_dsk = {
            (self._name, 0): (
                pd.Series,
                (
                    process_val_weights,
                    merged_key,
                    self.input_npartitions,
                    (self._name, 0, 0),
                ),
                qs,
                None,
                self.frame._meta.name,
            )
        }
        return {**dtype_dsk, **percentiles_dsk, **merge_dsk, **last_dsk}
