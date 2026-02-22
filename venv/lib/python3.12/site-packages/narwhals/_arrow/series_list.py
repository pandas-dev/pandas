from __future__ import annotations

from typing import TYPE_CHECKING

import pyarrow as pa
import pyarrow.compute as pc

from narwhals._arrow.utils import ArrowSeriesNamespace, list_agg, list_sort
from narwhals._compliant.any_namespace import ListNamespace
from narwhals._utils import not_implemented

if TYPE_CHECKING:
    from narwhals._arrow.series import ArrowSeries


class ArrowSeriesListNamespace(ArrowSeriesNamespace, ListNamespace["ArrowSeries"]):
    def len(self) -> ArrowSeries:
        return self.with_native(pc.list_value_length(self.native).cast(pa.uint32()))

    def get(self, index: int) -> ArrowSeries:
        return self.with_native(pc.list_element(self.native, index))

    def min(self) -> ArrowSeries:
        return self.with_native(list_agg(self.native, "min"))

    def max(self) -> ArrowSeries:
        return self.with_native(list_agg(self.native, "max"))

    def mean(self) -> ArrowSeries:
        return self.with_native(list_agg(self.native, "mean"))

    def median(self) -> ArrowSeries:
        return self.with_native(list_agg(self.native, "approximate_median"))

    def sum(self) -> ArrowSeries:
        return self.with_native(list_agg(self.native, "sum"))

    def sort(self, *, descending: bool, nulls_last: bool) -> ArrowSeries:
        return self.with_native(
            list_sort(self.native, descending=descending, nulls_last=nulls_last)
        )

    unique = not_implemented()
    contains = not_implemented()
