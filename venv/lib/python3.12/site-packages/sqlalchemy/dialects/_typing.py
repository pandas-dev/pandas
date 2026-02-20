# dialects/_typing.py
# Copyright (C) 2005-2026 the SQLAlchemy authors and contributors
# <see AUTHORS file>
#
# This module is part of SQLAlchemy and is released under
# the MIT License: https://www.opensource.org/licenses/mit-license.php
from __future__ import annotations

from typing import Any
from typing import Iterable
from typing import Mapping
from typing import Optional
from typing import Union

from ..sql import roles
from ..sql.base import ColumnCollection
from ..sql.schema import Column
from ..sql.schema import ColumnCollectionConstraint
from ..sql.schema import Index


_OnConflictConstraintT = Union[str, ColumnCollectionConstraint, Index, None]
_OnConflictIndexElementsT = Optional[
    Iterable[Union[Column[Any], str, roles.DDLConstraintColumnRole]]
]
_OnConflictIndexWhereT = Optional[roles.WhereHavingRole]
_OnConflictSetT = Optional[
    Union[Mapping[Any, Any], ColumnCollection[Any, Any]]
]
_OnConflictWhereT = Optional[roles.WhereHavingRole]
