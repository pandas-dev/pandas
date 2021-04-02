from __future__ import annotations

from collections import defaultdict

# from contextlib import contextmanager
# import copy
# from functools import partial
from typing import (
    Any,
    Callable,
    DefaultDict,
    Dict,
    List,
    Optional,
    Sequence,
    Tuple,
    Union,
)

# from uuid import uuid4
# import warnings
#
# import numpy as np
#
# from pandas._config import get_option
#
# from pandas._libs import lib
# from pandas._typing import (
#     Axis,
#     FrameOrSeries,
#     FrameOrSeriesUnion,
#     IndexLabel,
# )
# from pandas.compat._optional import import_optional_dependency
# from pandas.util._decorators import doc
#
# from pandas.core.dtypes.generic import ABCSeries
#
# import pandas as pd
# from pandas.api.types import is_list_like
# from pandas.core import generic
# import pandas.core.common as com
# from pandas.core.frame import DataFrame
# from pandas.core.generic import NDFrame
# from pandas.core.indexes.api import Index

BaseFormatter = Union[str, Callable]
ExtFormatter = Union[BaseFormatter, Dict[Any, Optional[BaseFormatter]]]
CSSPair = Tuple[str, Union[str, int, float]]
CSSList = List[CSSPair]
CSSProperties = Union[str, CSSList]
CSSStyles = List[Dict[str, CSSProperties]]  # = List[CSSDict]
# class CSSDict(TypedDict):  # available when TypedDict is valid in pandas
#     selector: str
#     props: CSSProperties


class StylerRender:
    def __init__(self):
        # variables will be overwritten by parent class
        self.hidden_index: bool = False
        self.hidden_columns: Sequence[int] = []
        self.ctx: DefaultDict[Tuple[int, int], CSSList] = defaultdict(list)
        self.cell_context: DefaultDict[Tuple[int, int], str] = defaultdict(str)
        self._todo: List[Tuple[Callable, Tuple, Dict]] = []

    def render(self, **kwargs) -> str:
        """
        Render the ``Styler`` including all applied styles to HTML.

        Parameters
        ----------
        **kwargs
            Any additional keyword arguments are passed
            through to ``self.template.render``.
            This is useful when you need to provide
            additional variables for a custom template.

        Returns
        -------
        rendered : str
            The rendered HTML.

        Notes
        -----
        Styler objects have defined the ``_repr_html_`` method
        which automatically calls ``self.render()`` when it's the
        last item in a Notebook cell. When calling ``Styler.render()``
        directly, wrap the result in ``IPython.display.HTML`` to view
        the rendered HTML in the notebook.

        Pandas uses the following keys in render. Arguments passed
        in ``**kwargs`` take precedence, so think carefully if you want
        to override them:

        * head
        * cellstyle
        * body
        * uuid
        * table_styles
        * caption
        * table_attributes
        """
        self._compute()
        # TODO: namespace all the pandas keys
        d = self._translate()
        d.update(kwargs)
        return self.template.render(**d)

    def _compute(self):
        """
        Execute the style functions built up in `self._todo`.

        Relies on the conventions that all style functions go through
        .apply or .applymap. The append styles to apply as tuples of

        (application method, *args, **kwargs)
        """
        self.ctx.clear()
        r = self
        for func, args, kwargs in self._todo:
            r = func(self)(*args, **kwargs)
        return r
